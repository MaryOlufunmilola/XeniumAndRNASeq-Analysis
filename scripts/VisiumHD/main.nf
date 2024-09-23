#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = "data/*.h5"  // Adjust according to your input files
params.output = "results/"

process PreprocessVisium {
    tag "preprocess"

    input:
    path input_file

    output:
    path "${params.output}/preprocessed_${input_file.baseName}.RData"

    script:
    """
    Rscript -e "
    library(Seurat)
    library(SeuratData)
    
    # Load data
    visium_data <- Read10X_h5('${input_file}')
    
    # Create Seurat object
    seurat_obj <- CreateSeuratObject(counts = visium_data)
    
    # Normalize the data
    seurat_obj <- NormalizeData(seurat_obj)
    
    # Find variable features
    seurat_obj <- FindVariableFeatures(seurat_obj)
    
    # Save the processed object
    save(seurat_obj, file='${params.output}/preprocessed_${input_file.baseName}.RData')
    "
    """
}

process AnalyzeSeurat {
    tag "analyze"

    input:
    path seurat_data

    output:
    path "${params.output}/analysis_${seurat_data.baseName}.RData"

    script:
    """
    Rscript -e "
    library(Seurat)
    
    # Load the preprocessed object
    load('${seurat_data}')
    
    # Run PCA
    seurat_obj <- RunPCA(seurat_obj)
    
    # Find clusters
    seurat_obj <- FindNeighbors(seurat_obj)
    seurat_obj <- FindClusters(seurat_obj)
    
    # Save the analysis results
    save(seurat_obj, file='${params.output}/analysis_${seurat_data.baseName}.RData')
    "
    """
}

workflow {
    preprocessed_files = Channel.fromPath(params.input)
    preprocessed_results = PreprocessVisium(preprocessed_files)
    AnalyzeSeurat(preprocessed_results)
}
