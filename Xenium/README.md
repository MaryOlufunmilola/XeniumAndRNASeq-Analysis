**Xenium spatial transcriptomics analysis project using R**

This repository contains scripts and workflows for analyzing spatial transcriptomics data from the **10x Genomics Xenium** platform using **R**. It supports preprocessing, visualization, spatial clustering, and downstream biological interpretation.

## Features

- Load and explore Xenium output using `Seurat`
- Visualize spatial gene expression and cell segmentation
- Perform clustering and marker gene identification
- Conduct spatial statistics and neighborhood analysis

## Setup

### 1. Clone the Repository

```bash
git clone https://github.com/yourusername/xenium-analysis-r.git
cd xenium-analysis-r
````

### 2. Install Required R Packages

You can install the required packages manually:

```r
install.packages(c(
  "Seurat", "SpatialExperiment", "ggplot2", "patchwork", 
  "dplyr", "BiocManager", "sp", "sf", "spdep"
))

# From Bioconductor
BiocManager::install("spatialLIBD")
BiocManager::install("SingleCellExperiment")
```

## Project Structure

```
xenium-analysis-r/
├── data/                # Raw and processed data (not included)
├── scripts/             # R scripts for loading, preprocessing, plotting
├── notebooks/           # R Markdown or Quarto files for analysis
├── results/             # Output plots and tables
├── renv.lock            # R environment lock file
└── README.md
```

## Usage

### Load and Preprocess Xenium Data

```r
source("scripts/load_xenium_data.R")
```

### Visualize Gene Expression

```r
rmarkdown::run("notebooks/visualize_genes.Rmd")
```

### Cluster and Annotate Cells

```r
rmarkdown::run("notebooks/clustering_and_annotation.Rmd")
```

### Spatial Analysis

* Neighborhood enrichment
* Ligand-receptor interaction (optional)

## Example Output

* Spatial plots of gene expression
* Clustering results overlaid on tissue
* Differentially expressed genes per cluster

## References

* [10x Genomics Xenium Documentation](https://www.10xgenomics.com/products/xenium)
* [Seurat](https://satijalab.org/seurat/)
* [SpatialExperiment](https://bioconductor.org/packages/SpatialExperiment/)
* [spdep](https://cran.r-project.org/web/packages/spdep/index.html)

## License

MIT License

## Contact

Questions or contributions welcome! Please open an issue or email [maryfunmisan@gmail.com](mailto:maryfunmisan@gmail.com)
