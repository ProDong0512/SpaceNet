# SpaceNet

# A xxxxxxxxx 介绍

这里插入workflow图

This repository describes how to analyse your data with SpaceNet.

Details on SpaceNet: 
* Installation
* Usage
* Relavants

See also the associated publication in **xxxx**: xxxxxxxxxx


---
## Quick start

### Installation of SpaceNet
    #installation of conda_env, recommended version python3.8.20, consistent with development env
    cd SpaceNet-release
    conda env create -n spacenet python=3.8.20 
    conda activate spacenet

    #installation of requirements
    pip install -r requirements.txt

### Accessibility to other softwares/database (Optional)
* [CellTrek](https://github.com/navinlabcode/CellTrek)
* [Cell2location](https://github.com/BayraktarLab/cell2location)
* [Tangram](https://github.com/broadinstitute/Tangram)


We recommend using 
    [this notebook](notebooks/PBMC10k_SCENIC-protocol-CLI.ipynb) 
    as a template for running an interactive analysis in Jupyter.
See the 
    [installation instructions](docs/installation.md)
    for information on setting up a kernel with pySCENIC and other required packages.

### Running the Nextflow pipeline on the example dataset

#### Requirements (Nextflow/containers)

The following tools are required to run the steps in this Nextflow pipeline:
* [Nextflow](https://www.nextflow.io/)
* A container system, either of:
    * [Docker](https://docs.docker.com/)
    * [Singularity](https://www.sylabs.io/singularity/)

The following container images will be pulled by nextflow as needed:
* Docker: [aertslab/pyscenic:latest](https://hub.docker.com/r/aertslab/pyscenic).
* Singularity: [aertslab/pySCENIC:latest](https://www.singularity-hub.org/collections/2033).
* [See also here.](https://github.com/aertslab/pySCENIC#docker-and-singularity-images)

#### Using the test profile

A quick test can be accomplished using the `test` profile, which automatically pulls the testing dataset (described in full below):

    nextflow run aertslab/SCENICprotocol \
        -profile docker,test

This small test dataset takes approximately 70s to run using 6 threads on a standard desktop computer.

#### Download testing dataset

Alternately, the same data can be run with a more verbose approach (this is more illustrative for how to substitute other data into the pipeline).
Download a minimum set of SCENIC database files for a human dataset (approximately 78 MB).

    mkdir example && cd example/
    # Transcription factors:
    wget https://raw.githubusercontent.com/aertslab/SCENICprotocol/master/example/test_TFs_tiny.txt
    # Motif to TF annotation database:
    wget https://raw.githubusercontent.com/aertslab/SCENICprotocol/master/example/motifs.tbl
    # Ranking databases:
    wget https://raw.githubusercontent.com/aertslab/SCENICprotocol/master/example/genome-ranking.feather
    # Finally, get a tiny sample expression matrix (loom format):
    cd scNiche-main

    conda env create -f scniche_dev.yaml -n scniche
  conda activate scniche
      wget https://raw.githubusercontent.com/aertslab/SCENICprotocol/master/example/expr_mat_tiny.loom


#### Running the example pipeline

Either Docker or Singularity images can be used by specifying the appropriate profile (`-profile docker` or `-profile singularity`).
Please note that for the tiny test dataset to run successfully, the default thresholds need to be lowered.

##### Using loom input

    nextflow run aertslab/SCENICprotocol \
        -profile docker \
        --loom_input expr_mat_tiny.loom \
        --loom_output pyscenic_integrated-output.loom \
        --TFs test_TFs_tiny.txt \
        --motifs motifs.tbl \
        --db *feather \
        --thr_min_genes 1

By default, this pipeline uses the container specified by the `--pyscenic_container` parameter.
This is currently set to `aertslab/pyscenic:0.9.19`, which uses a container with both pySCENIC and Scanpy `1.4.4.post1` installed.
A custom container can be used (e.g. one built on a local machine) by passing the name of this container to the `--pyscenic_container` parameter.

##### Expected output

The output of this pipeline is a loom-formatted file (by default: `output/pyscenic_integrated-output.loom`) containing:
* The original expression matrix
* The pySCENIC-specific results:
    * Regulons (TFs and their target genes)
    * AUCell matrix (cell enrichment scores for each regulon)
    * Dimensionality reduction embeddings based on the AUCell matrix (t-SNE, UMAP)
*  Results from the parallel best-practices analysis using highly variable genes:
    * Dimensionality reduction embeddings (t-SNE, UMAP)
    * Louvain clustering annotations

## General requirements for this workflow
* Python version 3.6 or greater
* Tested on various Unix/Linux distributions (Ubuntu 18.04, CentOS 7.6.1810, MacOS 10.14.5)

---

## Relavants and more information




