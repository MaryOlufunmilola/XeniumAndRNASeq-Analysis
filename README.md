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
The downstream differential expression analysis is [here](https://github.com/MaryOlufunmilola/Bioinformatics-Workflows/blob/master/RNA-Seq/scripts/RNA-Seq.Rmd).
