# Repository of gene features used for gene prioritization

This repository contains gene features (gene x feature matrices) used in the Finucane lab, primarily for gene prioritization methods.

## Raw data

Most raw data is in the form of (cell x gene) or (tissue x gene) count or TPM matrices. Experiments were performed either at the bulk or single cell level, and in some cases are single cells are merged before we derive features. Raw data will in many cases not be available on the public repository due to either total size or because it is currently unpublished and we do not have permissions to share. Please update this list everytime a new dataset is added.

- [human_multiple](https://www.gtexportal.org/home/) from [the GTEx consortium](http://dx.doi.org/10.1038/nature24277)
- [human_immune](https://preview.data.humancellatlas.org) from Aviv Regev and colleagues
- [human_heme](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74246) from [Corces _et al._](https://www.nature.com/articles/ng.3646)
- [human_pbmc](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc8k) from [Zheng _et al._](https://www.nature.com/articles/ncomms14049)
- [human_pancreas](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85241) from [Murano _et al._](https://www.sciencedirect.com/science/article/pii/S2405471216302927)
- human_gut from Alex Shalek, Jose Ordovas-Montanes, and colleagues
- [human_brain](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930) from [Lake _et al._](http://science.sciencemag.org/content/352/6293/1586)
- [human_brain2](http://www.brainspan.org) from [Miller _et al._](https://www.nature.com/articles/nature13185)
- [mouse_immume](https://www.immgen.org) from Christophe Benoist, Hideyuki Yoshida, Caleb Lareau, and colleagues
- [mouse_brain](http://dropviz.org) from [Saunders _et al._](https://www.biorxiv.org/content/early/2018/04/20/299081)
- mouse_aorta from Raj Gupta, Aditya Kalluri, and colleagues

## Feature types

The gerenal goal is to derive several types of gene expression-based features:

- Features that explain the most variance across cell populations (gene loadings on top 100 PCs, ...)
- Features that represent genes that define predefined cell populations or identified clusters (one vs. all, gene expression profiles in clusters, ...)
- Features that represent expression programs shared across cell types (co-expression gene modules, ...)

Features can be found [here](https://github.com/FinucaneLab/gene_features/tree/master/features).

## Getting setup 

Run [install.R](https://github.com/FinucaneLab/gene_features/tree/master/code/install.R) to install all necessary packages. Working code for each type of data is named similarly and lives in [code](https://github.com/FinucaneLab/gene_features/tree/master/code/).

## Analysis

First steps:
- Read in, QC, filter, scale, and normalize data (i.e. [features/mouse\_brain/variablegenes.pdf](https://github.com/FinucaneLab/gene_features/tree/master/features/mouse_brain/variablegenes.pdf))

- Perform PCA across all cells (where available) or meta-cells or tissues (i.e. [features/mouse\_brain/PCA.pdf](https://github.com/FinucaneLab/gene_features/tree/master/features/mouse_brain/PCA.pdf))

Current output features:
- Unweighted gene loadings (i.e. [features/mouse\_brain/u\_matrix.txt](https://github.com/FinucaneLab/gene_features/tree/master/features/mouse_brain/u_matrix.txt))
- Normalized expression matrices for predefined or identified clusters (i.e. [features/mouse\_brain/ave\_expr.txt](https://github.com/FinucaneLab/gene_features/tree/master/features/mouse_brain/ave_expr.txt))

Next steps:
- t-test one vs. all (slow but could run in parallel)?
- WCGNA?
- Verify clusters without predifined biological classes
- UMAP (rather than t-SNE) visualizations of clusters


