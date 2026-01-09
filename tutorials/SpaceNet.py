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

def cal_spot_radius(data):
    train_neighbors = NearestNeighbors(n_neighbors=2, metric='euclidean').fit(data)
    _, idx = train_neighbors.kneighbors(data)
    adj = train_neighbors.kneighbors_graph(data, mode='distance')
    spot_radius = adj.toarray().sum(axis=1).mean()/2
    return spot_radius


def construct_circlegraph(centers, radius):
    if (type(radius) == int) or (type(radius) == float):
        circles = [shapely.geometry.Point(x).buffer(radius) for x in centers]
        circlegraph = shapely.union_all(circles)
    else:
        assert len(centers) == len(radius), 'radius must be a fixed value, or an array the same length with centers'
        circles = []
        for i in range(len(centers)):
            circles.append(shapely.geometry.Point(centers[i]).buffer(radius[i]))
        circlegraph = shapely.union_all(circles)
    return circlegraph

def generate_search_range(adata, radius=None, disease_metric_key=str, cut_off=False, search_factor=1):
    if cut_off:
        adata_ds = adata[adata.obs[disease_metric_key] > cut_off, :]
    else:
        adata_ds = adata[adata.obs[disease_metric_key] != 0, :]
    adata_ds.obs[f'{disease_metric_key}_norm'] = (adata_ds.obs[disease_metric_key]/max(adata_ds.obs[disease_metric_key])).values
    centers = adata_ds.obsm['spatial']
    if radius:
        radius = radius
    else:
        spot_radius = cal_spot_radius(adata.obsm['spatial'])
        radius = adata_ds.obs[f'{disease_metric_key}_norm'].values*spot_radius*(search_factor-1) + spot_radius
    cg = construct_circlegraph(centers, radius)
    return cg, radius

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
    
    query_dict, query_failed = query_protein_name(protein_list, species=species, db_path=db_path)
    query_dict_t = {v: k for k, v in query_dict.items()}
    string_ids = list(query_dict.values())
    query_result = query_interactions(string_ids=string_ids, species=species, db_path=db_path, score_threshold=score_threshold)
    query_result['node_1'] = [query_dict_t[x] for x in query_result['protein1']]
    query_result['node_2'] = [query_dict_t[x] for x in query_result['protein2']]
    query_result = query_result[['node_1', 'protein1', 'node_2', 'protein2', 'combined_score']]
    
    return query_result, query_failed

from scipy.sparse import issparse
def cal_gene_enrichment(adata, disease_metric_key=str, cut_off=0, log=True, mode='wilcoxon'):    
    if log:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    if mode == 'wilcoxon':
        adata.var['gene_enrichment'] = 0
        adata.obs[disease_metric_key] = adata.obs[disease_metric_key].astype('category')
        sc.tl.rank_genes_groups(adata, groupby=disease_metric_key, method='wilcoxon')
        result = adata.uns['rank_genes_groups']
        df_geneenrich = pd.DataFrame({'gene': [result['names'][x][1] for x in range(len(result['names']))], 
                              'gene_enrichment': [result['scores'][x][1] for x in range(len(result['scores']))], 
                              'pvals': [result['pvals'][x][1] for x in range(len(result['pvals']))], 
                              'pvals_adj': [result['pvals_adj'][x][1] for x in range(len(result['pvals_adj']))]}).set_index('gene')
        df_geneenrich = df_geneenrich.loc[adata.var_names]
        adata.var[df_geneenrich.columns] = df_geneenrich
        return adata
    elif mode == 'pcc':
        adata.var['gene_enrichment'] = 0
        if cut_off:
            adata_ds = adata[adata.obs[disease_metric_key] != 0, :]
        else:
            adata_ds = adata[adata.obs[disease_metric_key] > cut_off, :]
        disease_values = adata_ds.obs[disease_metric_key].values
        for gene in adata_ds.var_names:
            exp = adata_ds[:,gene].X
            if issparse(exp):
                exp = exp.A.flatten()
            pcc = np.corrcoef(exp, disease_values)[0, 1]
            adata.var['gene_enrichment'].loc[gene] = pcc
        return adata

def cal_enriched_genes(adata, disease_metric_key=str, pval=None, pval_adj=None, top_g=None, bottom_g=None):
    if pval_adj:
        enriched_genes = adata[:, adata.var['pvals_adj'] < pval_adj].var.copy()
        enriched_genes = enriched_genes.sort_values(by='pvals_adj', ascending=True)
        return enriched_genes
    elif pval:
        enriched_genes = adata[:, adata.var['pvals'] < pval].var.copy()
        enriched_genes = enriched_genes.sort_values(by='pvals', ascending=True)
        return enriched_genes        
    elif (top_g and bottom_g):
        enriched_genes = adata.var.copy()
        enriched_genes_top = enriched_genes.sort_values(by='gene_enrichment', ascending=False).head(top_g)
        enriched_genes_bottom = enriched_genes.sort_values(by='gene_enrichment', ascending=True).head(bottom_g)
        enriched_genes = pd.concat([enriched_genes_top, enriched_genes_bottom]).sort_values(by='gene_enrichment', ascending=False)
        return enriched_genes
    elif top_g:
        enriched_genes = adata.var.copy()
        enriched_genes = enriched_genes.sort_values(by='gene_enrichment', ascending=False).head(top_g)
        return enriched_genes   
    elif bottom_g:
        enriched_genes = adata.var.copy()
        enriched_genes = enriched_genes.sort_values(by='gene_enrichment', ascending=True).head(bottom_g)
        return enriched_genes
    else:
        raise Exception("Please give a pval/pval_adj or the number of top/bottom genes")

        
def cal_genecoexp(adata, query_result):  # method: 'pearson', 'kendall', 'spearman'
    query_result['coexpef'] = 0
    gene_list = adata.var_names.tolist()
    for i in range(len(query_result)):
        gene_idx1 = gene_list.index(query_result.iloc[i]['node_1'])
        exp1 = adata.X[:, gene_idx1]
        gene_idx2 = gene_list.index(query_result.iloc[i]['node_2'])
        exp2 = adata.X[:, gene_idx2]
        if issparse(exp1):
            exp1 = exp1.A.flatten()
            exp2 = exp2.A.flatten()
        else:
            exp1 = exp1.flatten()
            exp2 = exp2.flatten()
        query_result.loc[i, 'coexpef'] =  np.corrcoef(exp1, exp2)[0, 1]
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
    for ct in list(expression_ratios.columns):
        for gene in list(expression_ratios.index):
            if selected_genes:
                if gene in selected_genes:
                    node1_list.append(gene)
                    node2_list.append(ct)
                    ratio_list.append(expression_ratios.loc[gene][ct])
                else:
                    continue
            else:
                node1_list.append(gene)
                node2_list.append(ct)
                ratio_list.append(expression_ratios.loc[gene][ct])
    gene_ct_net = pd.DataFrame({'node1':node1_list, 'node2':node2_list, 'ratio':ratio_list})
    return gene_ct_net