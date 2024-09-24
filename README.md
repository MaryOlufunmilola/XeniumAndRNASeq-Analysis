# RNA-Seq Data Analysis

Several bash commands are run on the Linux to download and manipulate RNA-Seq data. 

Using nfcore pipeline, **FastQC** was run to check the fastq data quality. Sequence reads were trimmed to remove possible adapter sequences and nucleotides with poor quality using **Trim Galore**. The trimmed reads were mapped to the Homo sapiens GRCh38 reference genome available on ENSEMBL using the **STAR** aligner. The STAR aligner is a splice aligner that detects splice junctions and incorporates them to help align the entire read sequences. Used **Salmon** to perform the downstream BAM-level quantification which calculates the number of unique reads per gene/exon.  

``` r
nextflow run \
    nf-core/rnaseq -r 3.15.0\
    --input "Samplesheet.csv" \
    --outdir "\/efs\/OUTDir\/" \
    --gtf "Homo_sapiens.GRCh38.112.gtf" \
    --fasta "Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa" \
    --skip_bigwig true \
    --skip_deseq2_qc true \
    -work-dir "\/efs\/work\/" \
    -profile docker
```

Using **DESeq2** package in R, a comparison of gene expression between treatment and control samples was performed. The Wald test was used to generate p-values and log2 fold changes. Genes with an adjusted p-value < 0.05 and absolute log2 fold change > 1 were called as differentially expressed genes for the comparison. 
**ClusterProfiler** was implemented to perform statistical analysis and visualization of functional profiles for genes and sets of genes.
The downstream differetial expression analysis is here.


# scRNA-Seq Data Analysis


# VisiumHD Data Analysis

#Nextflow pipeline for analyzing Visium HD data using the Seurat package in R
Nextflow, R and Seurat are the basic required software.
Nextflow Script named `main.nf` contains the lines of code to preprocess and analyze ((e.g., differential expression, visualization)) the Visium HD data using Seurat.
Visium HD `.h5` data files are in a directory named `data` and the preprocessed and analyzed data in the `results/` directory.

```bash
nextflow run main.nf -c nextflow.config
```


# Xenium Data Analysis
