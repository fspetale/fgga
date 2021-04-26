----------
# fgga

# R package

# FGGA: Factor Graph Gene ontology Annotation

FGGA is a graph-based machine learning approach for the automated and consistent GO annotation of protein coding genes. The input is a set of GO-term annotated protein coding genes previously characterized in terms of a fixed number of user-defined features, including the presence/absence of PFAM domains, physical-chemical properties, presence of signal peptides, among others. The set of GO-terms defines the output GO subgraph. A hierarchical ensemble (SVMs) machine learning model is generated. This model can be used to predict the GO subgraph annotations of uncharacterized protein coding genes. Individual GO-term annotations are accompanied by maximum a posteriori probability estimates issued by the native message passing algorithm of factor graphs.


# INSTALLATION

The fgga R source package can be directly downloaded from [Bioconductor repository](https://bioconductor.org/) or [GitHub repository](https://github.com/fspetale/fgga). This R package contains a experimental dataset as example, one pre-run R object and all functions needed to run FGGA.

\## From Bioconductor repository

if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")}
      
BiocManager::install("fgga")

\## Or from GitHub repository using devtools

BiocManager::install("devtools")

devtools::install_github("fspetale/fgga")

# REFERENCES

1: Spetale F.E., Tapia E., Krsticevic F., Roda F. and Bulacio P. “A Factor Graph Approach to Automated GO Annotation”. PLoS ONE 11(1): e0146986, 2016.

2: Spetale Flavio E., Arce D., Krsticevic F., Bulacio P. and Tapia E. “Consistent prediction of GO protein localization”. Scientific Report 7787(8), 2018
