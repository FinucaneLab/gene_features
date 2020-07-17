# Repository of gene features used for gene prioritization

This repository contains gene features (gene x feature matrices) used in the Finucane lab, primarily for gene prioritization methods.

## Raw data

Most raw data is in the form of (cell x gene) or (tissue x gene) count or TPM matrices. Experiments were performed either at the bulk or single cell level, and in some cases are single cells are merged before we derive features. Raw data will in many cases not be available on the public repository due to either total size or because it is currently unpublished and we do not have permissions to share. Please update this list everytime a new dataset is added.

- [human_multiple](https://www.gtexportal.org/home/) from [the GTEx consortium](https://www.biorxiv.org/content/10.1101/787903v1) (v8)
- [human_immune]((https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79)) from the Human Cell Atlas
- [human_heme](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74246) from [Corces _et al._](https://www.nature.com/articles/ng.3646)
- [human_pbmc](https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc8k) from [Zheng _et al._](https://www.nature.com/articles/ncomms14049)
- [human_pancreas](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE85241) from [Murano _et al._](https://www.sciencedirect.com/science/article/pii/S2405471216302927)
- [human_gut from](https://singlecell.broadinstitute.org/single_cell/study/SCP259/intra-and-inter-cellular-rewiring-of-the-human-colon-during-ulcerative-colitis) [Smilie _et _al.](https://doi.org/10.1016/j.cell.2019.06.029)
- [human_brain](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97930) from [Lake _et al._](http://science.sciencemag.org/content/352/6293/1586)
- [human_brain2](http://www.brainspan.org) from [Miller _et al._](https://www.nature.com/articles/nature13185)
- [mouse_immume](https://www.immgen.org) from [Yoshia _et al._](https://doi.org/10.1016/j.cell.2018.12.036)
- [mouse_brain](http://dropviz.org) from [Saunders _et al._](https://www.biorxiv.org/content/early/2018/04/20/299081)
- [mouse_aorta](https://singlecell.broadinstitute.org/single_cell/study/SCP289/single-cell-analysis-of-the-normal-mouse-aorta-reveals-functionally-distinct-endothelial-cell-populations) from [Kalluri _et al._](https://doi.org/10.1161/CIRCULATIONAHA.118.038362)
- [mouse_islets](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121416) from [Sharon _et al._](https://doi.org/10.1016/j.cell.2018.12.003)
- [mouse_digestive](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95630) from [Gao _et al._](https://www.nature.com/articles/s41556-018-0105-4)
- [mouse_heart](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118545) from [Hu _et al._](http://genesdev.cshlp.org/content/32/19-20/1344)
- [human_colon](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95459) from [Kinchen _et al._](https://doi.org/10.1016/j.cell.2018.08.067)
- [mouse_multiple](https://figshare.com/s/865e694ad06d5857db4b) from [Han _et al._](https://doi.org/10.1016/j.cell.2018.02.001)
- [mouse_lung](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119228) from [Cohen _et al._](https://doi.org/10.1016/j.cell.2018.09.009)
- [mouse_thymus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107910) from [Kernfeld _et al._](https://doi.org/10.1016/j.immuni.2018.04.015)
- [human_placenta](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89497) from [Liu _et al._](https://doi.org/10.1038/s41422-018-0066-y)
- [human_brain3](http://development.psychencode.org/files/processed_data/scRNA-seq/) from [Li _et al._](https://science.sciencemag.org/content/362/6420/eaat7615/)
- [mouse_microglia](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121654) from [Hammond _et al._](https://doi.org/10.1016/j.immuni.2018.11.004)
- [mouse_gutendoderm](https://endoderm-explorer.com/) from [Nowotschin _et al._](https://www.nature.com/articles/s41586-019-1127-1)
- [mouse_development](https://oncoscape.v3.sttrcancer.org/atlas.gs.washington.edu.mouse.rna/downloads) from [Cao _et al._](https://www.nature.com/articles/s41586-019-0969-x)
- [human_coloncancer](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146771) from [Zhang _et al._](https://doi.org/10.1016/j.cell.2020.03.048)
- [mouse_nerve](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142541) from [Wolbert _et al._](https://doi.org/10.1073/pnas.1912139117)
- [mouse_hemogenicendothelium](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137116) from [Zhu _et al._](https://doi.org/10.1182/blood.2020004801)
- [human_retina](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116106) from [Lu _et al._](https://doi.org/10.1016/j.devcel.2020.04.009)
- [human_pancreasductal](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131886) from [Qadir _et al._](https://doi.org/10.1073/pnas.1918314117)
- [mouse_hairfollicle](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115424) from [Shin _et al._](https://doi.org/10.1016/j.devcel.2020.03.019)
- [mouse_adipocyte](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145477) from [Zhong _et al._](https://elifesciences.org/articles/54695)
- [mouse_muscle](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143437) from [Micheli _et al._](https://doi.org/10.1016/j.celrep.2020.02.067)
- [mouse_airway](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8221) from [Miller _et al._](https://doi.org/10.1016/j.devcel.2020.01.033)
- [human_thymus](https://zenodo.org/record/3572422) from [Park _et al._](https://doi.org/10.1126/science.aay3224)
- [human_colon2](https://www.gutcellatlas.org/) from [James _et al._](https://doi.org/10.1038/s41590-020-0602-z)
- [mouse_endothelium](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8077) from [Kalucka _et al._](https://doi.org/10.1016/j.cell.2020.01.015)
- [mouse_vagina](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142212) from [Ali _et al._](https://doi.org/10.1016/j.celrep.2020.01.003)
- [human_eye](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135133) from [Orozco _et al._](10.1016/j.celrep.2019.12.082)
- [human_hippocampus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119212) from [Zhong _et al._](https://doi.org/10.1038/s41586-019-1917-5)
- [human_muscle](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130646) from [Rubenstein _et al._](https://doi.org/10.1038/s41598-019-57110-6)
- [human_intestine](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125970) from [Wang _et al._](https://doi.org/10.1084/jem.20191130)
- [human_kidney](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131685) from [Liao _et al._](https://doi.org/10.1038/s41597-019-0351-8)
- [human_monocytes](https://singlecell.broadinstitute.org/single_cell/study/SCP43/atlas-of-human-blood-dendritic-cells-and-monocytes) from [Villani _et al._](https://doi.org/10.1126/science.aah4573)
- [human_liver](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136103) from [Dobie _et al._](https://doi.org/10.1016/j.celrep.2019.10.024)
- [human_fetalblood](https://doi.org/10.1038/s41586-019-1652-y) from [Popescu _et al._](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-7407)
- [human_bonemarrow](https://doi.org/10.1172/jci.insight.124928) from [Oetjen _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120221)
- [mouse_brain2](https://doi.org/10.1126/science.aam8999) from [Rosenberg _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3017261)
~~- [mouse_brain3](https://doi.org/10.1016/j.cell.2018.06.021) from [Zeisel _et al._](https://http://mousebrain.org/downloads.html)~~
- [mouse_brain4](https://doi.org/10.1016/j.cell.2019.05.006) from [Welch _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126836)
- [human_brain4](https://doi.org/10.1016/j.cell.2019.05.006) from [Welch _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126836)
~~- [mouse_brain5](https://doi.org/10.1016/j.cell.2019.09.020) from [Kim _et al._](https://data.mendeley.com/datasets/ypx3sw2f7c/3)~~
- [mouse_kidney](https://doi.org/10.1016/j.devcel.2019.10.005) from [Rasnick _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129798)
- [human_nk](https://doi.org/10.1172/jci.insight.133103) from [Ferrari de Andrade _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139249)
- [human_retina2a](https://doi.org/10.1038/s41467-019-12780-8) from [Menon _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137537)
- [human_kidney2](https://doi.org/10.1126/science.aat5031) from [Stewart _et al._](https://data.humancellatlas.org/explore/projects/abe1a013-af7a-45ed-8c26-f3793c24a1f4)
- [human_tcell](https://doi.org/10.1038/s41467-019-12464-3) from [Szabo _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126030)
- [human_bladder](https://doi.org/10.1681/ASN.2019040335) from [Yu _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129845)
- [human_kidney3](https://doi.org/10.1073/pnas.1908706116) from [Wilson _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131882)
- [human_embryo](https://doi.org/10.1038/s41586-019-1500-0) from [Zhou _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109555)
- [mouse_epithelium](https://doi.org/10.1038/s41556-019-0378-2) from [Sharir _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131204)
- [human_lymphnodes](https://doi.org/10.1016/j.immuni.2019.06.027) from [Takeda _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124494)
- [mouse_multiple2](https://www.biorxiv.org/content/10.1101/661728v3) from [Pisco _et al._](https://figshare.com/articles/dataset/tms_gene_data/11413869)
- [human_lung](https://doi.org/10.1126/sciadv.aba1983) from [Adams _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831)
- [mouse_gastrulation](https://doi.org/10.1038/s41586-019-0933-9) from [Pijuan-Sala _et al._](https://content.cruk.cam.ac.uk/jmlab)
- [human_prostate](https://doi.org/10.1016/j.celrep.2018.11.086) from [Henry _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117403)
- [human_ileum](https://doi.org/10.1016/j.cell.2019.08.008) from [Martin _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134809)
- [human_airway](https://www.biorxiv.org/content/10.1101/2019.12.21.884759v1) from [Deprez _et al._](https://www.genomique.eu/cellbrowser/HCA/)
- [human_csf](https://doi.org/10.1038/s41467-019-14118-w) from [Schafflick _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138266)
- [human_testis](https://doi.org/10.1016/j.celrep.2018.10.026) from [Hermann _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109037)
- [human_synovialfibroblast](https://doi.org/10.1038/s41586-020-2222-z) from [Wei _et al._](https://singlecell.broadinstitute.org/single_cell/study/SCP469/)
- [mouse_muscle2](https://doi.org/10.1016/j.celrep.2020.02.067) from [De Micheli _et al._](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143435)

## Feature types

We derive features underlying tisse type / cellular processes using:

- Features that explain the most variance across (and within) cell populations (gene loadings on top X PCs, gene loading on top X ICs, within cluster PCs and ICs)
- Features that represent genes that define predefined cell populations or identified clusters (one vs. all differential gene expression -- t-stat, DE genes)
- Features that represent expression programs shared across cell types (TBD, co-expression gene modules)

Features can be found [here](https://github.com/FinucaneLab/gene_features/tree/master/features).

## Getting setup 

(Need to update) Run [install.R](https://github.com/FinucaneLab/gene_features/tree/master/code/install.R) to install all necessary packages. Working code for each type of data is named similarly and lives in [code](https://github.com/FinucaneLab/gene_features/tree/master/code/).

## Analysis

- Read in, QC, filter, scale, and normalize data (i.e. [plots/human\_immune/variablegenes.pdf](https://github.com/FinucaneLab/gene_features/tree/master/plots/human_immune/variablegenes.pdf))

- Perform PCA and ICA across all cells (where available) or meta-cells or tissues (i.e. [plots/human\_immune/pcaelbow.pdf](https://github.com/FinucaneLab/gene_features/tree/master/plots/human_immune/pcaelbow.pdf))

- Perform clustering and UMAP and plot features on projection (i.e. [plots/human\_immune/umap\_clusters.pdf](https://github.com/FinucaneLab/gene_features/tree/master/plots/human_immune/umap_clusters.pdf)

- Perform differential expression analysis (i.e. [plots/human\_immune/umap\_degenes.pdf](https://github.com/FinucaneLab/gene_features/tree/master/plots/human_immune/umap_degenes.pdf)

## Derived features

Current output features:
- Unweighted gene loadings from PCA (i.e. [features/human\_immune/projected\_pcaloadings.txt.gz](https://github.com/FinucaneLab/gene_features/tree/master/features/human_immune/projected_pcaloadings.txt.gz))
- Unweighted gene loadings from ICA (i.e. [features/human\_immune/projected\_icaloadings.txt.gz](https://github.com/FinucaneLab/gene_features/tree/master/features/human_immune/projected_icaloadings.txt.gz))
- Unweighted gene loadings from PCA within cluster (i.e. [features/human\_immune/projected\_pcaloadings\_clusters.txt.gz](https://github.com/FinucaneLab/gene_features/tree/master/features/human_immune/projected_pcaloadings\_clusters.txt.gz))
- Average expression across clusters (pre-defined and identified) (i.e. [features/human\_immune/average\_expression.txt](https://github.com/FinucaneLab/gene_features/tree/master/features/human_immune/average_expression.txt))
- One-vs-all t-test statistic (pre-defined and identified) (i.e. [features/human\_immune/diffexprs\_tstat\_clusters.txt](https://github.com/FinucaneLab/gene_features/tree/master/features/human_immune/diffexprs_tstat_clusters.txt))
- Differentially expressed genes across clusters (pre-defined and identified) (i.e. [features/human\_immune/diffexprs\_genes\_clusters.txt](https://github.com/FinucaneLab/gene_features/tree/master/features/human_immune/diffexprs_genes_clusters.txt))

Next steps:
- Gene modules (WCGNA)
- Co-expression analysis

