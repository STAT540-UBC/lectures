
# Clustering and Differential Analysis

## 1. Graph-based clustering

### A. KNN graph construction and community detection (PhenoGraph)

- **Levine et al.** (2015) [*Cell* 162:184-197](https://doi.org/10.1016/j.cell.2015.05.047). PhenoGraph: KNN graph + Jaccard-weighted SNN + Louvain community detection, applied to mass cytometry of AML.

### Further reading

- **Newman & Girvan** (2004) [*Phys Rev E* 69:026113](https://doi.org/10.1103/PhysRevE.69.026113). Finding and evaluating community structure in networks. Defines the modularity function Q.

### B. Modularity optimization: Louvain and Leiden

- **Blondel et al.** (2008) [*J Stat Mech* 2008:P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008). Fast unfolding of communities in large networks.
- **Traag et al.** (2019) [*Sci Rep* 9:5233](https://doi.org/10.1038/s41598-019-41695-z). From Louvain to Leiden: guaranteeing well-connected communities.

### C. Graph Laplacian and spectral clustering

- **Wang et al.** (2017) [*Nat Methods* 14:414-416](https://doi.org/10.1038/nmeth.4207). SIMLR: multi-kernel similarity learning + spectral clustering for single-cell RNA-seq visualization and analysis.

### D. Stochastic block models

- **Morelli et al.** (2021) [*BMC Bioinformatics* 22:576](https://doi.org/10.1186/s12859-021-04489-7). Nested stochastic block models applied to the analysis of single cell data. Principled alternative to modularity optimization with hierarchical structure and statistical model selection.


## 2. Challenges and significance of clustering

### A. Review of clustering challenges

- **Kiselev et al.** (2019) [*Nat Rev Genet* 20:273-282](https://www.nature.com/articles/s41576-018-0088-9). Challenges in unsupervised clustering of single-cell RNA-seq data.

### B. Significance analysis for clustering

- **Grabski et al.** (2023) [*Nat Methods* 20(8):1196-1202](https://doi.org/10.1038/s41592-023-01933-9). Significance analysis for clustering with single-cell RNA-sequencing data.

### C. Inference after latent variable estimation (double dipping)

- **Neufeld et al.** (2024) [*Biostatistics* 25(1):270-287](https://doi.org/10.1093/biostatistics/kxac047). Inference after latent variable estimation for single-cell RNA sequencing data. Proposes count splitting to avoid double dipping.

### D. Avoiding over-clustering

- **DenAdel et al.** (2025) [*Am J Hum Genet* 112(4):940-951](https://doi.org/10.1016/j.ajhg.2025.02.014). Artificial variables help to avoid over-clustering in single-cell RNA sequencing.
