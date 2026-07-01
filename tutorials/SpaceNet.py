import warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os
import random
import shapely
import duckdb
import time
from sklearn.neighbors import NearestNeighbors
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from scipy.sparse import issparse
from scipy.stats import pearsonr




def cal_spot_radius(data):
    """
    Estimate the physical radius of one spatial spot.

    The radius is defined as half of the average nearest-neighbor
    distance among spatial spots.
    """
    data = np.asarray(data)

    if data.shape[0] < 2:
        raise ValueError("At least two spatial spots are required to estimate spot radius.")

    train_neighbors = NearestNeighbors(
        n_neighbors=2,
        metric='euclidean'
    ).fit(data)

    distances, _ = train_neighbors.kneighbors(data)

    # distances[:, 0] is the distance to itself, usually 0.
    # distances[:, 1] is the nearest-neighbor distance.
    spot_radius = np.mean(distances[:, 1]) / 2

    return float(spot_radius)



def construct_circlegraph(centers, radius):
    """
    Construct a union polygon from circular spot footprints.

    Parameters
    ----------
    centers : array-like, shape (n_spots, 2)
        Spatial coordinates of selected spots.
    radius : float or array-like
        If float, all spots use the same radius.
        If array-like, each spot uses its own radius.

    Returns
    -------
    circlegraph : shapely geometry
        Union of all spot circles.
    """
    centers = np.asarray(centers)

    if len(centers) == 0:
        return shapely.geometry.GeometryCollection()

    if np.isscalar(radius):
        circles = [
            shapely.geometry.Point(x).buffer(float(radius))
            for x in centers
        ]
    else:
        radius = np.asarray(radius)

        assert len(centers) == len(radius), (
            "radius must be a fixed value, or an array with the same length as centers"
        )

        circles = [
            shapely.geometry.Point(centers[i]).buffer(float(radius[i]))
            for i in range(len(centers))
        ]

    circlegraph = shapely.union_all(circles)

    return circlegraph



def generate_search_range(
    adata,
    radius=None,
    disease_metric_key=None,
    cut_off=None,
    search_factor=None,
    conservative=True,
    return_mask=False
):
    """
    Generate pathology-associated spatial search range.

    Revised conservative default:
    - Select pathology-positive spots using disease_metric_key.
    - If cut_off is None, spots with disease_metric != 0 are selected.
    - If cut_off is provided, spots with disease_metric > cut_off are selected.
    - By default, no signal-dependent expansion is applied.
      Each selected spot contributes only its physical spot footprint.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics object with adata.obsm['spatial'].
    radius : float or array-like or None
        Fixed radius for spot footprints. If None, spot radius is estimated
        from spatial coordinates.
    disease_metric_key : str
        Column in adata.obs used to define pathology-positive spots.
    cut_off : float or None
        Optional threshold for filtering background signal.
        If None, disease_metric != 0 is used.
        If float, disease_metric > cut_off is used.
    search_factor : float or None
        Deprecated for main analysis. Only used when conservative=False.
        In conservative=True mode, search_factor is ignored.
    conservative : bool
        If True, use fixed spot radius without signal-dependent expansion.
        If False, use the original signal-weighted expansion behavior.
    return_mask : bool
        If True, also return roi_mask for selected spots.

    Returns
    -------
    cg : shapely geometry
        Union of selected spot footprints.
    radius_out : float or numpy.ndarray
        Radius used for ROI construction.
    roi_mask : numpy.ndarray, optional
        Boolean mask indicating selected ROI seed spots.
    """
    if disease_metric_key is None:
        raise ValueError("Please provide disease_metric_key.")

    if disease_metric_key not in adata.obs.columns:
        raise ValueError(f"{disease_metric_key} not found in adata.obs.")

    disease_values = adata.obs[disease_metric_key].astype(float).values

    # ROI seed selection
    if cut_off is None:
        # No artificial cutoff: use nonzero pathology-positive spots.
        roi_mask = disease_values != 0
    else:
        # Optional background filtering.
        roi_mask = disease_values > cut_off

    centers = adata.obsm['spatial'][roi_mask]

    if len(centers) == 0:
        cg = shapely.geometry.GeometryCollection()
        if return_mask:
            return cg, np.nan, roi_mask
        return cg, np.nan

    spot_radius = cal_spot_radius(adata.obsm['spatial'])

    # Radius definition
    if radius is not None:
        radius_out = radius

    elif conservative:
        # Revised main-analysis behavior:
        # no signal-dependent radius expansion.
        radius_out = spot_radius

    else:
        # Original signal-weighted expansion behavior.
        # Retained only for sensitivity analysis or backward compatibility.
        if search_factor is None:
            search_factor = 1

        selected_values = disease_values[roi_mask]
        max_value = np.nanmax(selected_values)

        if max_value <= 0 or np.isnan(max_value):
            disease_norm = np.ones(len(selected_values))
        else:
            disease_norm = selected_values / max_value

        radius_out = disease_norm * spot_radius * (search_factor - 1) + spot_radius

    cg = construct_circlegraph(centers, radius_out)

    if return_mask:
        return cg, radius_out, roi_mask

    return cg, radius_out



def _align_c2l_to_adata(
    c2l_df,
    adata,
    celltype_cols=None
):
    """
    Align Cell2location spot-by-cell-type matrix to adata.obs_names.

    Parameters
    ----------
    c2l_df : pandas.DataFrame
        Cell2location result matrix. Rows are spatial spots and columns are cell types.
    adata : AnnData
        Spatial transcriptomics object.
    celltype_cols : list or None
        Cell-type columns to use. If None, all numeric columns are used.

    Returns
    -------
    c2l_aligned : pandas.DataFrame
        Cell2location matrix aligned to adata.obs_names.
    celltype_cols : list
        Cell-type columns used for enrichment analysis.
    """
    c2l_df = c2l_df.copy()

    common_spots = adata.obs_names.intersection(c2l_df.index)

    if len(common_spots) == 0:
        raise ValueError(
            "No overlapping spot names between adata.obs_names and c2l_df.index. "
            "Please check barcode format."
        )

    if len(common_spots) < adata.n_obs:
        print(
            f"Warning: only {len(common_spots)} / {adata.n_obs} adata spots "
            "are found in c2l_df. Subsetting adata-aligned C2L matrix."
        )

    c2l_aligned = c2l_df.loc[adata.obs_names.intersection(c2l_df.index)].copy()

    if celltype_cols is None:
        celltype_cols = c2l_aligned.select_dtypes(include=[np.number]).columns.tolist()
    else:
        missing_cols = [x for x in celltype_cols if x not in c2l_aligned.columns]
        if len(missing_cols) > 0:
            raise ValueError(f"These celltype_cols are missing from c2l_df: {missing_cols}")

    c2l_aligned[celltype_cols] = c2l_aligned[celltype_cols].astype(float)

    return c2l_aligned, celltype_cols



def cal_celltype_enrichment_c2l(
    c2l_df,
    adata,
    roi_mask,
    celltype_cols=None,
    anatomy_key=None,
    n_permutations=1000,
    exclude_roi_in_null=True,
    min_total_abundance=1e-9,
    random_state=0
):
    """
    Calculate ROI-associated cell-type enrichment from Cell2location output.

    This function is the Cell2location-compatible version of cal_celltype_enrichment().
    Instead of counting mapped cells, it sums spot-level cell-type abundance/proportion.

    Parameters
    ----------
    c2l_df : pandas.DataFrame
        Cell2location result matrix. Rows are spatial spots and columns are cell types.
    adata : AnnData
        Spatial transcriptomics object. Requires adata.obs_names and optionally anatomy_key.
    roi_mask : array-like of bool
        Boolean vector indicating ROI spots. Length must equal adata.n_obs.
    celltype_cols : list or None
        Cell-type columns in c2l_df. If None, all numeric columns are used.
    anatomy_key : str or None
        Optional adata.obs column for anatomy-matched background sampling.
    n_permutations : int
        Number of random background samplings.
    exclude_roi_in_null : bool
        Whether to exclude true ROI spots from background sampling.
    min_total_abundance : float
        Cell types with total abundance below this threshold are flagged.
    random_state : int
        Random seed.

    Returns
    -------
    result_df : pandas.DataFrame
        Cell-type enrichment table.
    c2l_aligned : pandas.DataFrame
        C2L matrix aligned to the adata spots used in the analysis.
    """
    rng = np.random.default_rng(random_state)

    roi_mask = np.asarray(roi_mask).astype(bool)

    if len(roi_mask) != adata.n_obs:
        raise ValueError("roi_mask length must equal adata.n_obs")

    c2l_aligned, celltype_cols = _align_c2l_to_adata(
        c2l_df=c2l_df,
        adata=adata,
        celltype_cols=celltype_cols
    )

    # If some adata spots are missing from C2L, subset roi_mask and anatomy labels accordingly.
    aligned_spots = c2l_aligned.index
    adata_spot_index = pd.Index(adata.obs_names)
    aligned_indices = adata_spot_index.get_indexer(aligned_spots)

    roi_mask_aligned = roi_mask[aligned_indices]

    X = c2l_aligned[celltype_cols].values.astype(float)

    anatomy_labels = None
    if anatomy_key is not None:
        if anatomy_key not in adata.obs.columns:
            raise ValueError(f"{anatomy_key} not found in adata.obs")
        anatomy_labels = adata.obs.iloc[aligned_indices][anatomy_key].astype(str).values

    # Observed abundance
    observed_abundance_roi = X[roi_mask_aligned, :].sum(axis=0)
    total_abundance = X.sum(axis=0)

    total_abundance_all_celltypes = float(X.sum())
    roi_abundance_all_celltypes = float(X[roi_mask_aligned, :].sum())

    # Permutation null abundance
    null_abundance = np.zeros((n_permutations, len(celltype_cols)), dtype=float)

    for i in range(n_permutations):
        sampled_spots = _sample_matched_background_spots(
            roi_mask=roi_mask_aligned,
            anatomy_labels=anatomy_labels,
            n_spots=int(roi_mask_aligned.sum()),
            exclude_roi=exclude_roi_in_null,
            rng=rng
        )

        null_abundance[i, :] = X[sampled_spots, :].sum(axis=0)

    rows = []

    for j, ct in enumerate(celltype_cols):
        observed = float(observed_abundance_roi[j])
        total = float(total_abundance[j])

        null = null_abundance[:, j]
        expected = float(np.mean(null))
        null_sd = float(np.std(null, ddof=1))

        observed_fraction = observed / total if total > 0 else np.nan
        expected_fraction = expected / total if total > 0 else np.nan

        fold_enrichment = (observed + 1e-9) / (expected + 1e-9)

        p_enrich = (np.sum(null >= observed) + 1) / (n_permutations + 1)
        p_deplete = (np.sum(null <= observed) + 1) / (n_permutations + 1)
        p_two_sided = min(1.0, 2 * min(p_enrich, p_deplete))

        rows.append({
            "cell_type": ct,
            "abundance_roi": observed,
            "abundance_total": total,
            "total_abundance_all_celltypes": total_abundance_all_celltypes,
            "roi_abundance_all_celltypes": roi_abundance_all_celltypes,
            "n_roi_spots": int(roi_mask_aligned.sum()),
            "n_total_spots": int(len(roi_mask_aligned)),
            "roi_spot_fraction": float(roi_mask_aligned.mean()),
            "observed_fraction": observed_fraction,
            "expected_abundance_roi": expected,
            "expected_fraction": expected_fraction,
            "null_sd": null_sd,
            "fold_enrichment": fold_enrichment,
            "p_enrich": p_enrich,
            "p_deplete": p_deplete,
            "p_two_sided": p_two_sided,
            "low_abundance_flag": total < min_total_abundance
        })

    result_df = pd.DataFrame(rows)

    result_df["padj_enrich"] = multipletests(
        result_df["p_enrich"].values,
        method="fdr_bh"
    )[1]

    result_df["padj_two_sided"] = multipletests(
        result_df["p_two_sided"].values,
        method="fdr_bh"
    )[1]

    result_df = result_df.sort_values(
        by=["padj_enrich", "fold_enrichment"],
        ascending=[True, False]
    )

    return result_df, c2l_aligned




def cal_celltype_enrichment_c2l_multi(
    samples,
    celltype_cols=None,
    anatomy_key=None,
    n_permutations=1000,
    exclude_roi_in_null=True,
    min_total_abundance=1e-9,
    random_state=0
):
    """
    Pooled Cell2location cell-type enrichment across multiple spatial samples.

    Parameters
    ----------
    samples : list of dict
        Each dict should contain:
        {
            "sample": sample name,
            "adata": AnnData,
            "c2l_df": Cell2location DataFrame,
            "roi_mask": boolean ROI mask
        }
    celltype_cols : list or None
        Cell types to analyze. If None, use the union of numeric columns across samples.
    anatomy_key : str or None
        Optional anatomy-matched background sampling key.
    n_permutations : int
        Number of pooled permutations.
    exclude_roi_in_null : bool
        Whether to exclude true ROI spots from background sampling.
    min_total_abundance : float
        Low-abundance flag threshold.
    random_state : int
        Random seed.

    Returns
    -------
    pooled_result : pandas.DataFrame
        Pooled enrichment result across samples.
    sample_detail : pandas.DataFrame
        Observed ROI and total abundance per sample and cell type.
    """
    rng = np.random.default_rng(random_state)

    # Determine celltype columns
    if celltype_cols is None:
        ct_set = set()
        for s in samples:
            numeric_cols = s["c2l_df"].select_dtypes(include=[np.number]).columns.tolist()
            ct_set.update(numeric_cols)
        celltype_cols = sorted(ct_set)

    prepared = []
    sample_rows = []

    for s in samples:
        sample_name = s.get("sample", "sample")
        adata = s["adata"]
        c2l_df = s["c2l_df"]
        roi_mask = np.asarray(s["roi_mask"]).astype(bool)

        if len(roi_mask) != adata.n_obs:
            raise ValueError(f"roi_mask length mismatch in {sample_name}")

        c2l_aligned = c2l_df.copy()

        common_spots = adata.obs_names.intersection(c2l_aligned.index)
        if len(common_spots) == 0:
            raise ValueError(f"No overlapping spots in {sample_name}")

        c2l_aligned = c2l_aligned.loc[common_spots].copy()

        # Fill missing cell types with 0 so all samples share the same columns
        for ct in celltype_cols:
            if ct not in c2l_aligned.columns:
                c2l_aligned[ct] = 0.0

        c2l_aligned = c2l_aligned[celltype_cols].astype(float)

        adata_spot_index = pd.Index(adata.obs_names)
        aligned_indices = adata_spot_index.get_indexer(common_spots)
        roi_mask_aligned = roi_mask[aligned_indices]

        anatomy_labels = None
        if anatomy_key is not None:
            if anatomy_key not in adata.obs.columns:
                raise ValueError(f"{anatomy_key} not found in {sample_name}.obs")
            anatomy_labels = adata.obs.iloc[aligned_indices][anatomy_key].astype(str).values

        X = c2l_aligned.values.astype(float)

        observed = X[roi_mask_aligned, :].sum(axis=0)
        total = X.sum(axis=0)

        for j, ct in enumerate(celltype_cols):
            sample_rows.append({
                "sample": sample_name,
                "cell_type": ct,
                "abundance_roi": float(observed[j]),
                "abundance_total": float(total[j]),
                "n_roi_spots": int(roi_mask_aligned.sum()),
                "n_total_spots": int(len(roi_mask_aligned))
            })

        prepared.append({
            "sample": sample_name,
            "X": X,
            "roi_mask": roi_mask_aligned,
            "anatomy_labels": anatomy_labels,
            "observed": observed,
            "total": total
        })

    observed_pooled = np.sum([x["observed"] for x in prepared], axis=0)
    total_pooled = np.sum([x["total"] for x in prepared], axis=0)

    null_pooled = np.zeros((n_permutations, len(celltype_cols)), dtype=float)

    for i in range(n_permutations):
        pooled_iter = np.zeros(len(celltype_cols), dtype=float)

        for x in prepared:
            sampled_spots = _sample_matched_background_spots(
                roi_mask=x["roi_mask"],
                anatomy_labels=x["anatomy_labels"],
                n_spots=int(x["roi_mask"].sum()),
                exclude_roi=exclude_roi_in_null,
                rng=rng
            )

            pooled_iter += x["X"][sampled_spots, :].sum(axis=0)

        null_pooled[i, :] = pooled_iter

    rows = []

    for j, ct in enumerate(celltype_cols):
        observed = float(observed_pooled[j])
        total = float(total_pooled[j])

        null = null_pooled[:, j]
        expected = float(np.mean(null))
        null_sd = float(np.std(null, ddof=1))

        observed_fraction = observed / total if total > 0 else np.nan
        expected_fraction = expected / total if total > 0 else np.nan

        fold_enrichment = (observed + 1e-9) / (expected + 1e-9)

        p_enrich = (np.sum(null >= observed) + 1) / (n_permutations + 1)
        p_deplete = (np.sum(null <= observed) + 1) / (n_permutations + 1)
        p_two_sided = min(1.0, 2 * min(p_enrich, p_deplete))

        rows.append({
            "cell_type": ct,
            "abundance_roi": observed,
            "abundance_total": total,
            "observed_fraction": observed_fraction,
            "expected_abundance_roi": expected,
            "expected_fraction": expected_fraction,
            "null_sd": null_sd,
            "fold_enrichment": fold_enrichment,
            "p_enrich": p_enrich,
            "p_deplete": p_deplete,
            "p_two_sided": p_two_sided,
            "n_samples": len(samples),
            "low_abundance_flag": total < min_total_abundance
        })

    pooled_result = pd.DataFrame(rows)

    pooled_result["padj_enrich"] = multipletests(
        pooled_result["p_enrich"].values,
        method="fdr_bh"
    )[1]

    pooled_result["padj_two_sided"] = multipletests(
        pooled_result["p_two_sided"].values,
        method="fdr_bh"
    )[1]

    pooled_result = pooled_result.sort_values(
        by=["padj_enrich", "fold_enrichment"],
        ascending=[True, False]
    )

    sample_detail = pd.DataFrame(sample_rows)

    return pooled_result, sample_detail




def _sample_matched_background_spots(
    roi_mask,
    anatomy_labels=None,
    n_spots=None,
    exclude_roi=True,
    rng=None
):
    """
    Sample background spots matched to ROI spot number and optionally anatomy composition.
    """
    if rng is None:
        rng = np.random.default_rng()

    roi_mask = np.asarray(roi_mask).astype(bool)
    all_indices = np.arange(len(roi_mask))

    if n_spots is None:
        n_spots = int(roi_mask.sum())

    if n_spots <= 0:
        return np.array([], dtype=int)

    if anatomy_labels is None:
        if exclude_roi:
            pool = all_indices[~roi_mask]
        else:
            pool = all_indices

        if len(pool) == 0:
            raise ValueError("No background spots available for permutation sampling.")

        replace = len(pool) < n_spots
        return rng.choice(pool, size=n_spots, replace=replace)

    anatomy_labels = np.asarray(anatomy_labels)
    sampled = []

    for label in pd.Series(anatomy_labels[roi_mask]).value_counts().index:
        n_label_roi = int(np.sum((anatomy_labels == label) & roi_mask))

        if exclude_roi:
            pool = all_indices[(anatomy_labels == label) & (~roi_mask)]
        else:
            pool = all_indices[anatomy_labels == label]

        if len(pool) == 0:
            continue

        replace = len(pool) < n_label_roi
        sampled_label = rng.choice(pool, size=n_label_roi, replace=replace)
        sampled.extend(sampled_label.tolist())

    return np.array(sampled, dtype=int)



def _assign_cells_to_nearest_spot(
    cell_df,
    spot_coords,
    spot_names=None,
    x_col="coord_x",
    y_col="coord_y"
):
    """
    Assign each mapped cell to its nearest spatial transcriptomics spot.

    Parameters
    ----------
    cell_df : pandas.DataFrame
        CellTrek or mapped-cell table containing cell coordinates.
    spot_coords : array-like, shape (n_spots, 2)
        Spatial coordinates from adata.obsm['spatial'].
    spot_names : array-like, optional
        Spot identifiers, usually adata.obs_names.
    x_col, y_col : str
        Coordinate columns in cell_df.

    Returns
    -------
    cell_df_out : pandas.DataFrame
        Input cell table with nearest_spot_idx and nearest_spot added.
    """
    cell_df_out = cell_df.copy()

    if spot_names is None:
        spot_names = np.arange(spot_coords.shape[0])

    nbrs = NearestNeighbors(n_neighbors=1, metric="euclidean")
    nbrs.fit(spot_coords)

    cell_coords = cell_df_out[[x_col, y_col]].values
    _, indices = nbrs.kneighbors(cell_coords)

    nearest_idx = indices[:, 0]
    cell_df_out["nearest_spot_idx"] = nearest_idx
    cell_df_out["nearest_spot"] = np.array(spot_names)[nearest_idx]

    return cell_df_out





def query_protein_name(protein_list, species=None, db_path='/string/'):
    
    assert species in [None, 'Homo sapiens', 'Mus musculus', 'Rattus norvegicus'], 'Enter a valid species name(Homo sapiens, Mus musculus, Rattus norvegicus), or skip if unsure.'
    
    query_dict = {}
    query_failed = []
    
    if species == None:
        query_db = db_path + 'protein_db_aliases.duckdb'
    elif species == 'Homo sapiens':
        query_db = db_path + 'protein_aliases_human.duckdb'
    elif species == 'Mus musculus':
        query_db = db_path + 'protein_aliases_mouse.duckdb'
    elif species == 'Rattus norvegicus':
        query_db = db_path + 'protein_aliases_rat.duckdb'
    
    conn = duckdb.connect(query_db)
    id_list = ', '.join([f"'{id}'" for id in protein_list])
    
    query = f"""
        SELECT string_protein_id, alias
        FROM aliases
        WHERE alias IN ({id_list})
    """

    result = conn.sql(query).df()
    conn.close()
    result = result.drop_duplicates().set_index('alias')
    
    for protein in protein_list:
        try:
            if type(result.loc[protein]['string_protein_id']) == str:
                query_dict[protein] = result.loc[protein]['string_protein_id']
            else:
                query_dict[protein] = sorted(result.loc[protein]['string_protein_id'])[0]
            #print(query_dict[protein])
        except:
            query_failed.append(protein)
    
    return query_dict, query_failed



def query_interactions(string_ids: list, species=None, score_threshold=700, db_path='/string/'):
    
    assert species in [None, 'Homo sapiens', 'Mus musculus', 'Rattus norvegicus'], 'Enter a valid species name(Homo sapiens, Mus musculus, Rattus norvegicus), or skip if unsure.'
    
    if not string_ids:
        return pd.DataFrame()
    
    if species == None:
        query_db = db_path + 'protein_db_links.duckdb'
    elif species == 'Homo sapiens':
        query_db = db_path + 'protein_links_human.duckdb'
    elif species == 'Mus musculus':
        query_db = db_path + 'protein_links_mouse.duckdb'
    elif species == 'Rattus norvegicus':
        query_db = db_path + 'protein_links_rat.duckdb'
    
    conn = duckdb.connect(query_db)
    conn.execute("PRAGMA threads=1")
    
    id_list = ', '.join([f"'{id}'" for id in string_ids])
    query = f"""
        SELECT protein1, protein2, combined_score
        FROM interactions
        WHERE combined_score >= {score_threshold}
          AND (protein1 IN ({id_list}) OR protein2 IN ({id_list}))
    """
    result = conn.sql(query).df()
    conn.close()

    target_set = set(string_ids)
    filtered = result[
        result['protein1'].isin(target_set) & 
        result['protein2'].isin(target_set)
    ]
    
    return filtered



def query(protein_list: list, species=None, score_threshold=700, db_path='/string/'):
    """
    Query STRING interactions among an input gene/protein list.

    Parameters
    ----------
    protein_list : list
        Gene symbols or protein aliases.
    species : str or None
        One of None, 'Homo sapiens', 'Mus musculus', 'Rattus norvegicus'.
    score_threshold : int
        STRING combined score cutoff. 700 corresponds to high confidence.
    db_path : str
        Path to local STRING duckdb files.

    Returns
    -------
    query_result : pandas.DataFrame
        Edge table with node_1, protein1, node_2, protein2, combined_score.
    query_failed : list
        Input genes that failed STRING ID mapping.
    """
    empty_cols = ['node_1', 'protein1', 'node_2', 'protein2', 'combined_score']

    query_dict, query_failed = query_protein_name(
        protein_list,
        species=species,
        db_path=db_path
    )

    if len(query_dict) == 0:
        return pd.DataFrame(columns=empty_cols), query_failed

    query_dict_t = {v: k for k, v in query_dict.items()}
    string_ids = list(query_dict.values())

    query_result = query_interactions(
        string_ids=string_ids,
        species=species,
        db_path=db_path,
        score_threshold=score_threshold
    )

    if query_result is None or len(query_result) == 0:
        return pd.DataFrame(columns=empty_cols), query_failed

    query_result = query_result.copy()
    query_result['node_1'] = [query_dict_t[x] for x in query_result['protein1']]
    query_result['node_2'] = [query_dict_t[x] for x in query_result['protein2']]
    query_result = query_result[empty_cols]

    return query_result, query_failed



def cal_gene_enrichment(
    adata,
    disease_metric_key=str,
    cut_off=0,
    log=True,
    mode='wilcoxon'
):
    
    """
    For mode='wilcoxon':
        Use this mode for binary ROI labels, such as 0/1 pathology regions.

    For mode='pcc':
        Use this mode for continuous pathology signals, such as plaque burden
        or normalized amyloid intensity.

    cut_off behavior in mode='pcc':
        - cut_off=None: use spots with disease_metric != 0.
        - cut_off=float: use spots with disease_metric > cut_off.

    The output columns are written to adata.var:
        - gene_enrichment
        - pvals
        - pvals_adj
    """
    if log:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    if mode == 'wilcoxon':
        adata.var['gene_enrichment'] = 0
        adata.obs[disease_metric_key] = adata.obs[disease_metric_key].astype('category')

        sc.tl.rank_genes_groups(
            adata,
            groupby=disease_metric_key,
            method='wilcoxon'
        )

        result = adata.uns['rank_genes_groups']

        df_geneenrich = pd.DataFrame({
            'gene': [result['names'][x][1] for x in range(len(result['names']))],
            'gene_enrichment': [result['scores'][x][1] for x in range(len(result['scores']))],
            'pvals': [result['pvals'][x][1] for x in range(len(result['pvals']))],
            'pvals_adj': [result['pvals_adj'][x][1] for x in range(len(result['pvals_adj']))]
        }).set_index('gene')

        df_geneenrich = df_geneenrich.loc[adata.var_names]
        adata.var[df_geneenrich.columns] = df_geneenrich

        return adata

    elif mode == 'pcc':
        # Initialize output columns for all genes
        adata.var['gene_enrichment'] = np.nan
        adata.var['pvals'] = np.nan
        adata.var['pvals_adj'] = np.nan

        # Select spots used for Pearson correlation
        if cut_off is None:
            spot_mask = adata.obs[disease_metric_key].values != 0
        else:
            spot_mask = adata.obs[disease_metric_key].values > cut_off

        adata_ds = adata[spot_mask, :].copy()

        disease_values = adata_ds.obs[disease_metric_key].astype(float).values

        r_values = []
        p_values = []
        genes = []

        for gene in adata_ds.var_names:
            exp = adata_ds[:, gene].X

            if issparse(exp):
                exp = exp.A.flatten()
            else:
                exp = np.asarray(exp).flatten()

            # Guard against constant expression or constant disease signal
            if np.std(exp) == 0 or np.std(disease_values) == 0:
                r = np.nan
                p = 1.0
            else:
                r, p = pearsonr(exp, disease_values)

                if np.isnan(r):
                    r = np.nan
                    p = 1.0

            genes.append(gene)
            r_values.append(r)
            p_values.append(p)

        p_values_for_fdr = np.array(p_values, dtype=float)
        p_values_for_fdr[np.isnan(p_values_for_fdr)] = 1.0

        pvals_adj = multipletests(
            p_values_for_fdr,
            method='fdr_bh'
        )[1]

        df_geneenrich = pd.DataFrame({
            'gene': genes,
            'gene_enrichment': r_values,
            'pvals': p_values,
            'pvals_adj': pvals_adj
        }).set_index('gene')

        adata.var.loc[df_geneenrich.index, 'gene_enrichment'] = df_geneenrich['gene_enrichment']
        adata.var.loc[df_geneenrich.index, 'pvals'] = df_geneenrich['pvals']
        adata.var.loc[df_geneenrich.index, 'pvals_adj'] = df_geneenrich['pvals_adj']

        return adata

    else:
        raise ValueError("mode must be either 'wilcoxon' or 'pcc'")
        


        
def cal_enriched_genes(
    adata,
    disease_metric_key=str,
    pval=None,
    pval_adj=None,
    top_g=None,
    bottom_g=None,
    sort_by="gene_enrichment"
):
    """
    Select ROI/pathology-associated genes.

    Revised behavior:
    - First apply optional p-value or FDR filtering.
    - Then optionally select top and/or bottom genes by effect size.
    """
    enriched_genes = adata.var.copy()

    if pval_adj is not None:
        if "pvals_adj" not in enriched_genes.columns:
            raise ValueError("adata.var does not contain 'pvals_adj'. Run cal_gene_enrichment first.")
        enriched_genes = enriched_genes[enriched_genes["pvals_adj"] < pval_adj]

    elif pval is not None:
        if "pvals" not in enriched_genes.columns:
            raise ValueError("adata.var does not contain 'pvals'. Run cal_gene_enrichment first.")
        enriched_genes = enriched_genes[enriched_genes["pvals"] < pval]

    if len(enriched_genes) == 0:
        return enriched_genes

    if sort_by not in enriched_genes.columns:
        raise ValueError(f"{sort_by} not found in adata.var.")

    if top_g is not None and bottom_g is not None:
        enriched_genes_top = enriched_genes.sort_values(
            by=sort_by,
            ascending=False
        ).head(top_g)

        enriched_genes_bottom = enriched_genes.sort_values(
            by=sort_by,
            ascending=True
        ).head(bottom_g)

        enriched_genes = pd.concat([
            enriched_genes_top,
            enriched_genes_bottom
        ]).drop_duplicates()

        enriched_genes = enriched_genes.sort_values(
            by=sort_by,
            ascending=False
        )

        return enriched_genes

    elif top_g is not None:
        return enriched_genes.sort_values(
            by=sort_by,
            ascending=False
        ).head(top_g)

    elif bottom_g is not None:
        return enriched_genes.sort_values(
            by=sort_by,
            ascending=True
        ).head(bottom_g)

    else:
        if pval_adj is not None:
            return enriched_genes.sort_values(by="pvals_adj", ascending=True)
        elif pval is not None:
            return enriched_genes.sort_values(by="pvals", ascending=True)
        else:
            raise Exception("Please give a pval/pval_adj or the number of top/bottom genes.")
        
        
        

        
def _get_gene_expression_vector(adata, gene, spot_mask=None):
    """
    Extract one gene expression vector from AnnData.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics object.
    gene : str
        Gene name.
    spot_mask : array-like of bool or None
        Optional spot mask. If provided, correlation is calculated only on selected spots.

    Returns
    -------
    exp : numpy.ndarray or None
        Gene expression vector. None if gene is not found.
    """
    if gene not in adata.var_names:
        return None

    gene_idx = adata.var_names.get_loc(gene)
    exp = adata.X[:, gene_idx]

    if issparse(exp):
        exp = exp.A.flatten()
    else:
        exp = np.asarray(exp).flatten()

    if spot_mask is not None:
        exp = exp[np.asarray(spot_mask).astype(bool)]

    return exp


def cal_genecoexp(
    adata,
    query_result,
    spot_mask=None,
    adjust_method='fdr_bh',
    min_spots=3
):
    """
    Calculate spatial co-expression statistics for STRING-supported gene pairs.

    For each STRING edge, this function computes the Pearson correlation between
    the two genes' expression profiles across spatial transcriptomic spots.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics object.
    query_result : pandas.DataFrame
        STRING edge table returned by query(). Must contain node_1 and node_2.
    spot_mask : array-like of bool or None
        Optional mask to restrict co-expression calculation to selected spots,
        such as ROI spots. If None, all spots are used.
    adjust_method : str
        Multiple-testing correction method used by statsmodels.multipletests.
    min_spots : int
        Minimum number of spots required for Pearson correlation.

    Returns
    -------
    query_result : pandas.DataFrame
        Edge table with added columns:
        coexpef, coexp_pvalue, coexp_padj.
    """
    query_result = query_result.copy()

    if len(query_result) == 0:
        query_result['coexpef'] = []
        query_result['coexp_pvalue'] = []
        query_result['coexp_padj'] = []
        return query_result

    if 'node_1' not in query_result.columns or 'node_2' not in query_result.columns:
        raise ValueError("query_result must contain 'node_1' and 'node_2' columns")

    if spot_mask is not None:
        spot_mask = np.asarray(spot_mask).astype(bool)
        if len(spot_mask) != adata.n_obs:
            raise ValueError("spot_mask length must equal adata.n_obs")

    r_values = []
    p_values = []

    for i in range(len(query_result)):
        gene1 = query_result.iloc[i]['node_1']
        gene2 = query_result.iloc[i]['node_2']

        exp1 = _get_gene_expression_vector(adata, gene1, spot_mask=spot_mask)
        exp2 = _get_gene_expression_vector(adata, gene2, spot_mask=spot_mask)

        if exp1 is None or exp2 is None:
            r = np.nan
            p = 1.0
        elif len(exp1) < min_spots or len(exp2) < min_spots:
            r = np.nan
            p = 1.0
        elif np.std(exp1) == 0 or np.std(exp2) == 0:
            r = np.nan
            p = 1.0
        else:
            r, p = pearsonr(exp1, exp2)

            if np.isnan(r) or np.isnan(p):
                r = np.nan
                p = 1.0

        r_values.append(r)
        p_values.append(p)

    p_values_for_fdr = np.asarray(p_values, dtype=float)
    p_values_for_fdr[np.isnan(p_values_for_fdr)] = 1.0

    p_adj = multipletests(
        p_values_for_fdr,
        method=adjust_method
    )[1]

    query_result['coexpef'] = r_values
    query_result['coexp_pvalue'] = p_values
    query_result['coexp_padj'] = p_adj

    return query_result

def calculate_gene_expression_by_celltype(adata, celltype_column='cell_type'):
    
    if celltype_column not in adata.obs.columns:
        raise ValueError(f"'{celltype_column}' not in adata.obs")
    
    if issparse(adata.X):
        expression_matrix = adata.X.toarray()
    else:
        expression_matrix = adata.X
    
    cell_types = adata.obs[celltype_column].values
    unique_cell_types = np.unique(cell_types)
    gene_names = adata.var_names
    
    expression_sum_df = pd.DataFrame(index=gene_names, columns=unique_cell_types)
    expression_ratio_df = pd.DataFrame(index=gene_names, columns=unique_cell_types)
    
    for cell_type in unique_cell_types:
        cell_mask = cell_types == cell_type
        cell_subset = expression_matrix[cell_mask, :]
        
        gene_sums = np.sum(cell_subset, axis=0)
        expression_sum_df[cell_type] = gene_sums
    
    total_expression_per_gene = expression_sum_df.sum(axis=1)
    
    for cell_type in unique_cell_types:
        expression_ratio_df[cell_type] = expression_sum_df[cell_type] / total_expression_per_gene
    
    expression_ratio_df = expression_ratio_df.fillna(0)
    
    return expression_sum_df, expression_ratio_df

def construct_gene_ct_net(expression_ratios, selected_genes=None):
    node1_list = []
    node2_list = []
    ratio_list = []

    if selected_genes is not None:
        selected_genes = set(selected_genes)

    for ct in list(expression_ratios.columns):
        for gene in list(expression_ratios.index):
            if selected_genes is None or gene in selected_genes:
                node1_list.append(gene)
                node2_list.append(ct)
                ratio_list.append(expression_ratios.loc[gene, ct])

    gene_ct_net = pd.DataFrame({
        'node1': node1_list,
        'node2': node2_list,
        'ratio': ratio_list
    })

    return gene_ct_net


def _largest_connected_component_size(edge_df, nodes=None):
    """
    Calculate largest connected component size from an undirected edge table.
    """
    if nodes is None:
        nodes = set()
    else:
        nodes = set(nodes)

    adjacency = {}

    for _, row in edge_df.iterrows():
        a = row['node_1']
        b = row['node_2']

        nodes.add(a)
        nodes.add(b)

        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)

    if len(nodes) == 0:
        return 0

    visited = set()
    largest_size = 0

    for node in nodes:
        if node in visited:
            continue

        stack = [node]
        component = set()

        while stack:
            current = stack.pop()
            if current in component:
                continue

            component.add(current)

            for neighbor in adjacency.get(current, set()):
                if neighbor not in component:
                    stack.append(neighbor)

        visited.update(component)
        largest_size = max(largest_size, len(component))

    return largest_size



def summarize_ppi_network(
    edge_df,
    input_genes=None,
    coexp_col='coexpef',
    coexp_padj_col='coexp_padj',
    coexp_padj_cutoff=0.05
):
    """
    Summarize STRING-prior spatial co-expression network.

    Parameters
    ----------
    edge_df : pandas.DataFrame
        Edge table returned by query() and cal_genecoexp().
    input_genes : list or None
        Original input gene set. If provided, density is calculated using this size.
    coexp_col : str
        Column containing spatial co-expression coefficient.
    coexp_padj_col : str
        Column containing FDR-adjusted co-expression P values.
    coexp_padj_cutoff : float
        Significance cutoff for spatially supported edges.

    Returns
    -------
    metrics : dict
        Network-level metrics.
    """
    if edge_df is None:
        edge_df = pd.DataFrame()

    edge_df = edge_df.copy()

    if input_genes is not None:
        input_genes = list(dict.fromkeys(input_genes))
        n_input_genes = len(input_genes)
    else:
        if len(edge_df) == 0:
            input_genes = []
        else:
            input_genes = sorted(set(edge_df['node_1']).union(set(edge_df['node_2'])))
        n_input_genes = len(input_genes)

    max_possible_edges = n_input_genes * (n_input_genes - 1) / 2
    n_string_edges = len(edge_df)

    if len(edge_df) > 0 and coexp_padj_col in edge_df.columns:
        sig_edge_df = edge_df[edge_df[coexp_padj_col] < coexp_padj_cutoff].copy()
    else:
        sig_edge_df = pd.DataFrame(columns=edge_df.columns)

    n_spatial_sig_edges = len(sig_edge_df)

    if len(edge_df) > 0 and coexp_col in edge_df.columns:
        mean_abs_coexp_all = float(np.nanmean(np.abs(edge_df[coexp_col].values)))
        mean_coexp_all = float(np.nanmean(edge_df[coexp_col].values))
    else:
        mean_abs_coexp_all = np.nan
        mean_coexp_all = np.nan

    if len(sig_edge_df) > 0 and coexp_col in sig_edge_df.columns:
        mean_abs_coexp_sig = float(np.nanmean(np.abs(sig_edge_df[coexp_col].values)))
        mean_coexp_sig = float(np.nanmean(sig_edge_df[coexp_col].values))
    else:
        mean_abs_coexp_sig = np.nan
        mean_coexp_sig = np.nan

    if len(edge_df) > 0 and 'combined_score' in edge_df.columns:
        mean_string_score = float(np.nanmean(edge_df['combined_score'].values))
    else:
        mean_string_score = np.nan

    metrics = {
        'n_input_genes': n_input_genes,
        'n_string_edges': n_string_edges,
        'string_edge_density': (
            n_string_edges / max_possible_edges
            if max_possible_edges > 0 else np.nan
        ),
        'n_spatial_sig_edges': n_spatial_sig_edges,
        'spatial_sig_edge_density': (
            n_spatial_sig_edges / max_possible_edges
            if max_possible_edges > 0 else np.nan
        ),
        'mean_string_score': mean_string_score,
        'mean_coexp_all_edges': mean_coexp_all,
        'mean_abs_coexp_all_edges': mean_abs_coexp_all,
        'mean_coexp_sig_edges': mean_coexp_sig,
        'mean_abs_coexp_sig_edges': mean_abs_coexp_sig,
        'largest_component_size_all_edges': _largest_connected_component_size(
            edge_df,
            nodes=input_genes
        ),
        'largest_component_size_sig_edges': _largest_connected_component_size(
            sig_edge_df,
            nodes=input_genes
        )
    }

    return metrics


def _get_detected_genes(adata, min_cells=1):
    """
    Return genes detected in at least min_cells spots.
    """
    X = adata.X

    if issparse(X):
        detected_counts = np.asarray((X > 0).sum(axis=0)).flatten()
    else:
        detected_counts = np.sum(np.asarray(X) > 0, axis=0)

    genes = np.asarray(adata.var_names)
    return genes[detected_counts >= min_cells].tolist()


def _get_mean_expression(adata):
    """
    Calculate mean expression of each gene across spots.
    """
    X = adata.X

    if issparse(X):
        mean_expr = np.asarray(X.mean(axis=0)).flatten()
    else:
        mean_expr = np.mean(np.asarray(X), axis=0)

    return pd.Series(mean_expr, index=adata.var_names)


def _sample_expression_matched_genes(
    observed_genes,
    background_genes,
    mean_expr,
    n_bins=10,
    rng=None
):
    """
    Sample random genes matched to observed genes by average expression bins.

    The expression bins are constructed using both observed and background genes,
    but random genes are sampled only from the background gene set.
    """
    if rng is None:
        rng = np.random.default_rng()

    observed_genes = list(dict.fromkeys([
        g for g in observed_genes
        if g in mean_expr.index
    ]))

    background_genes = list(dict.fromkeys([
        g for g in background_genes
        if g in mean_expr.index
    ]))

    if len(observed_genes) == 0:
        raise ValueError("No observed genes available for expression-matched sampling.")

    if len(background_genes) == 0:
        raise ValueError("No background genes available for expression-matched sampling.")

    # Build bins using observed + background genes
    bin_genes = list(dict.fromkeys(observed_genes + background_genes))

    bin_df = pd.DataFrame({
        "gene": bin_genes,
        "mean_expr": mean_expr.loc[bin_genes].values
    }).dropna()

    # Robust quantile binning
    if bin_df["mean_expr"].nunique() <= 1:
        bin_df["expr_bin"] = 0
    else:
        expr_rank = bin_df["mean_expr"].rank(method="first")
        bin_df["expr_bin"] = pd.qcut(
            expr_rank,
            q=min(n_bins, len(bin_df)),
            labels=False,
            duplicates="drop"
        )

    gene_to_bin = bin_df.set_index("gene")["expr_bin"].to_dict()

    background_df = bin_df[bin_df["gene"].isin(background_genes)].copy()

    selected = []
    selected_set = set()

    for gene in observed_genes:
        if gene not in gene_to_bin:
            continue

        bin_id = gene_to_bin[gene]

        pool = background_df[
            (background_df["expr_bin"] == bin_id) &
            (~background_df["gene"].isin(selected_set))
        ]["gene"].values

        # If the matched bin is empty, fall back to all unused background genes
        if len(pool) == 0:
            pool = background_df[
                ~background_df["gene"].isin(selected_set)
            ]["gene"].values

        if len(pool) == 0:
            break

        sampled_gene = rng.choice(pool)
        selected.append(sampled_gene)
        selected_set.add(sampled_gene)

    # Fill if needed
    if len(selected) < len(observed_genes):
        remaining_pool = [
            g for g in background_genes
            if g not in selected_set
        ]

        n_needed = len(observed_genes) - len(selected)

        if len(remaining_pool) >= n_needed:
            extra = rng.choice(
                remaining_pool,
                size=n_needed,
                replace=False
            ).tolist()
        else:
            extra = rng.choice(
                background_genes,
                size=n_needed,
                replace=True
            ).tolist()

        selected.extend(extra)

    return selected[:len(observed_genes)]




def cal_ppi_network_null_model(
    adata,
    observed_genes,
    species=None,
    db_path='/string/',
    score_threshold=700,
    spot_mask=None,
    background_genes=None,
    min_cells=1,
    n_permutations=1000,
    expression_matched=True,
    n_expression_bins=10,
    coexp_padj_cutoff=0.05,
    random_state=0
):
    """
    Compare an observed STRING-prior spatial co-expression network against
    random gene-set networks.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics object.
    observed_genes : list
        ROI/pathology-associated genes used to build the observed network.
    species : str or None
        STRING species.
    db_path : str
        Path to local STRING duckdb files.
    score_threshold : int
        STRING combined score cutoff.
    spot_mask : array-like of bool or None
        Optional spot mask for spatial co-expression calculation.
    background_genes : list or None
        Background gene universe. If None, detected genes in adata are used.
    min_cells : int
        Minimum spots with nonzero expression for background gene inclusion.
    n_permutations : int
        Number of random networks.
    expression_matched : bool
        Whether to sample random genes matched by mean expression.
    n_expression_bins : int
        Number of expression bins for matched sampling.
    coexp_padj_cutoff : float
        FDR cutoff for spatially significant edges.
    random_state : int
        Random seed.

    Returns
    -------
    observed_edges : pandas.DataFrame
        Observed STRING-prior spatial co-expression edge table.
    observed_metrics : pandas.DataFrame
        One-row table of observed network metrics.
    null_metrics : pandas.DataFrame
        Metrics from random networks.
    null_pvalues : pandas.DataFrame
        Empirical P values and FDR across tested network metrics.
    """
    rng = np.random.default_rng(random_state)

    observed_genes = list(dict.fromkeys([
        g for g in observed_genes
        if g in adata.var_names
    ]))

    if len(observed_genes) < 2:
        raise ValueError("observed_genes must contain at least two genes found in adata.var_names")

    if background_genes is None:
        background_genes = _get_detected_genes(adata, min_cells=min_cells)
    else:
        background_genes = [
            g for g in background_genes
            if g in adata.var_names
        ]

    # Remove observed genes from random background when possible
    random_background = [
        g for g in background_genes
        if g not in set(observed_genes)
    ]

    if len(random_background) < len(observed_genes):
        random_background = background_genes

    mean_expr = _get_mean_expression(adata)

    # Observed network
    observed_edges, observed_query_failed = query(
        protein_list=observed_genes,
        species=species,
        score_threshold=score_threshold,
        db_path=db_path
    )

    observed_edges = cal_genecoexp(
        adata,
        observed_edges,
        spot_mask=spot_mask
    )

    observed_metrics_dict = summarize_ppi_network(
        observed_edges,
        input_genes=observed_genes,
        coexp_padj_cutoff=coexp_padj_cutoff
    )

    observed_metrics = pd.DataFrame([observed_metrics_dict])
    observed_metrics['n_query_failed'] = len(observed_query_failed)

    # Null networks
    null_rows = []

    for i in range(n_permutations):
        if expression_matched:
            random_genes = _sample_expression_matched_genes(
                observed_genes=observed_genes,
                background_genes=random_background,
                mean_expr=mean_expr,
                n_bins=n_expression_bins,
                rng=rng
            )
        else:
            replace = len(random_background) < len(observed_genes)
            random_genes = rng.choice(
                random_background,
                size=len(observed_genes),
                replace=replace
            ).tolist()

        random_edges, random_query_failed = query(
            protein_list=random_genes,
            species=species,
            score_threshold=score_threshold,
            db_path=db_path
        )

        random_edges = cal_genecoexp(
            adata,
            random_edges,
            spot_mask=spot_mask
        )

        random_metrics = summarize_ppi_network(
            random_edges,
            input_genes=random_genes,
            coexp_padj_cutoff=coexp_padj_cutoff
        )

        random_metrics['iteration'] = i
        random_metrics['n_query_failed'] = len(random_query_failed)

        null_rows.append(random_metrics)

    null_metrics = pd.DataFrame(null_rows)

    # Empirical P values for metrics where larger values indicate stronger network structure
    test_metrics = [
        'n_string_edges',
        'string_edge_density',
        'n_spatial_sig_edges',
        'spatial_sig_edge_density',
        'mean_abs_coexp_all_edges',
        'mean_abs_coexp_sig_edges',
        'largest_component_size_all_edges',
        'largest_component_size_sig_edges'
    ]

    pvalue_rows = []

    for metric in test_metrics:
        observed_value = observed_metrics_dict.get(metric, np.nan)

        if metric not in null_metrics.columns:
            continue

        null_values = null_metrics[metric].values.astype(float)
        null_values = null_values[~np.isnan(null_values)]

        if np.isnan(observed_value) or len(null_values) == 0:
            empirical_p = np.nan
        else:
            empirical_p = (
                np.sum(null_values >= observed_value) + 1
            ) / (
                len(null_values) + 1
            )

        pvalue_rows.append({
            'metric': metric,
            'observed': observed_value,
            'null_mean': float(np.nanmean(null_values)) if len(null_values) > 0 else np.nan,
            'null_sd': float(np.nanstd(null_values, ddof=1)) if len(null_values) > 1 else np.nan,
            'empirical_pvalue': empirical_p
        })

    null_pvalues = pd.DataFrame(pvalue_rows)

    if len(null_pvalues) > 0:
        p_for_fdr = null_pvalues['empirical_pvalue'].fillna(1).values
        null_pvalues['empirical_padj'] = multipletests(
            p_for_fdr,
            method='fdr_bh'
        )[1]

    return observed_edges, observed_metrics, null_metrics, null_pvalues
