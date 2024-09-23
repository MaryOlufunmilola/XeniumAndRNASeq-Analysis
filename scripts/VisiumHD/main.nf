---
title: "R Notebook"
output: html_notebook
---

```{r}
setwd("H:/Backup_I_J_drive/ScottNess/ACC")
library(AnnotationHub)
library(data.table)
library(scDblFinder)
library(qs)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(future)
library(arrow)
library(hdf5r)
library(tzdb)
library(BPCells)
library(Polychrome)
library(scCustomize)

n <-50
set.seed(72)
mycolors = as.vector(createPalette(n,  c('#FF0000', "#0000FF", '#FFFF00'), target = "normal", M=20000))
save(mycolors, file = "savefolder/colors.RData")
pdf("figurefolder/GeneralColorSet.pdf", height = 9, width = 13)
swatch(mycolors)
dev.off()

# needs to be set for large dataset analysis
options(Seurat.object.assay.version = "v5")
no_cores <- availableCores() - 4 ##setting options for future

if (!file.exists("figurefolder")){ dir.create("figurefolder") }  
if (!file.exists("filefolder")){ dir.create("filefolder") }  
if (!file.exists("savefolder")){ dir.create("savefolder") }  

localdir <- "HDvisium/"
files = list.files(localdir)
seurat_obj = list()
for (mfile in files){
    nfile <- gsub("sample_","",mfile)
    ss <- Load10X_Spatial(data.dir = paste0(localdir, mfile, "/binned_outputs/square_008um"), assay = "Spatial", slice = nfile)
    ss$orig.ident <- nfile
    ss$log10GenesPerUMI <- log10(ss$nFeature_Spatial) / log10(ss$nCount_Spatial)
    ss$spikeRatio <- PercentageFeatureSet(object = ss, pattern = "^ERCC")
	  ss$mitoRatio <- PercentageFeatureSet(object = ss, pattern = "^MT-")
    ss$MYBRatio <- PercentageFeatureSet(object = ss, pattern = "MYB$")
    ss <- NormalizeData(ss)
    ss <- FindVariableFeatures(ss)
    ss <- ScaleData(ss)
    saveRDS(ss, paste0("savefolder/",nfile,"_seu.rds"))
    
    pdf(paste0("figurefolder/VlnPlotinitialSpatial_", nfile,".pdf"), height = 9, width = 13)
    print(VlnPlot(ss, features = c("nFeature_Spatial", "nCount_Spatial", "mitoRatio", "MYBRatio"), ncol=2, cols = "lightsteelblue3", layer = "counts"))
    dev.off()
    vln.plot <-  VlnPlot(ss, features = "nCount_Spatial", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
    count.plot <- SpatialFeaturePlot(ss, features = "nCount_Spatial") + theme(legend.position = "right")
    pdf(paste0("figurefolder/VlnCountPlotinitialSpatial_", nfile,".pdf"), height = 9, width = 13)
    print(vln.plot | count.plot)
    dev.off()
    
    # we select 20,000 cells and create a new 'sketch' assay
    ss <- SketchData(ss,ncells = 20000, method = "LeverageScore",sketched.assay = "sketch")
    # switch analysis to sketched cells
    DefaultAssay(ss) <- "sketch"
    ss <- FindVariableFeatures(ss)
    ss <- ScaleData(ss)
    ss <- RunPCA(ss, assay = "sketch", reduction.name = "PCAsketch")
    ss <- FindNeighbors(ss, assay = "sketch", reduction = "PCAsketch", dims = 1:20)
    plan(multisession, workers = no_cores)
    options(future.globals.maxSize = 230 * 1024^3, future.rng.onMisuse = 'ignore') #default settings
    ss <- FindClusters(ss, cluster.name = "seurat_cluster_sketched", resolution = 0.5)
    ss <- RunUMAP(ss, reduction = "PCAsketch", reduction.name = "UMAPsketch", return.model = T, dims = 1:20)
    saveRDS(ss, paste0("savefolder/",nfile,"_seuske.rds"))
    
    Idents(ss) <- "seurat_cluster_sketched"
    p1 <- DimPlot(ss, reduction = "UMAPsketch", label = F) + ggtitle("Sketched clustering (20,000 cells)") + theme(legend.position = "bottom")
    
    ss <- ProjectData(ss, assay = "Spatial", full.reduction = "fullPCAsketch", sketched.assay = "sketch", sketched.reduction = "PCAsketch", umap.model = "UMAPsketch", dims = 1:20, refdata = list(seurat_cluster_projected = "seurat_cluster_sketched"))    
    saveRDS(ss, paste0("savefolder/",nfile,"_seufinal.rds"))
    
    # switch to full dataset
    DefaultAssay(ss) <- "Spatial"
    Idents(ss) <- "seurat_cluster_projected"
    p2 <- DimPlot(ss, reduction = "full.UMAPsketch", label = F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")
    
    pdf(paste0("figurefolder/UMAP_sketch_full_clusters_", nfile,".pdf"), height = 9, width = 13)
    print(p1 | p2)
    dev.off()
    
    pdf(paste0("figurefolder/UMAP_full_clusters_", nfile,".pdf"), height = 9, width = 13)
    print(DimPlot_scCustom(ss, reduction = "full.UMAPsketch", label = F, colors_use=mycolors, raster=FALSE) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom"))
    dev.off()
    
    pdf(paste0("figurefolder/Spatial_Allclusters_", nfile,".pdf"), height = 9, width = 13)
    print(SpatialDimPlot(ss, label = T, repel = T, label.size = 4))
    dev.off()
    
    #Highlight the spatial localization of a few clusters 
    #cells <- CellsByIdentities(ss, idents = levels(ss))
    #pdf(paste0("figurefolder/Spatial_Eachcluster_", nfile,".pdf"), height = 9, width = 13)
    #print(SpatialDimPlot(ss,cells.highlight = cells[setdiff(names(cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend())
    #dev.off()
    
    seurat_obj[[nfile]] <- ss
    #Find and visualize the top gene expression markers for each cluster:
    # Create down sampled object to make visualization either
    object_subset <- subset(ss, cells = Cells(ss[["Spatial"]]), downsample = 1000)
    
    # Order clusters by similarity
    DefaultAssay(object_subset) <- "Spatial"
    Idents(object_subset) <- "seurat_cluster_projected"
    object_subset <- BuildClusterTree(object_subset, assay = "Spatial", reduction = "fullPCAsketch", reorder = T)
    
    markers <- FindAllMarkers(object_subset, assay = "Spatial", only.pos = TRUE)
    write.csv(markers, file = paste0("filefolder/PositiveMarkers_downsampled", nfile,".pdf"))
    markers %>%
      group_by(cluster) %>%
      dplyr::filter(avg_log2FC > 1) %>%
      slice_head(n = 5) %>%
      ungroup() -> top5
    
    object_subset <- ScaleData(object_subset, assay = "Spatial", features = top5$gene)
    
    pdf(paste0("figurefolder/Heatmap_subset_", nfile,".pdf"), height = 9, width = 13)
    print(DoHeatmap(object_subset, assay = "Spatial", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend())
    dev.off()
}
rm(ss)
gc()
saveRDS(seurat_obj, "savefolder/Allsamples_Normdata.rds")

```

Analyze data across multiple samples.

```{r}
Vf_2 <- c(VariableFeatures(seurat_obj[[1]]), VariableFeatures(seurat_obj[[3]]))
Vf_5 <- c(VariableFeatures(seurat_obj[[2]]), VariableFeatures(seurat_obj[[4]]))
seurat_obj[["less<2"]] <- merge(x = seurat_obj[[1]], y = seurat_obj[[3]], add.cell.id = c('ACC251463','ACC368571'), merge.data =TRUE)
seurat_obj[["more>5"]] <- merge(x = seurat_obj[[2]], y = seurat_obj[[4]], add.cell.id = c('ACC317286','ACC742434'), merge.data =TRUE)
VariableFeatures(seurat_obj[["less<2"]]) <- Vf_2
VariableFeatures(seurat_obj[["more>5"]]) <- Vf_5

seurat_obj[1:4] <- NULL
saveRDS(seurat_obj[["more>5"]], "savefolder/samplesmorethan5years.rds")
saveRDS(seurat_obj[["less<2"]], "savefolder/samplesLessthan2years.rds")

ss <- seurat_obj[[nfile]]
ss <- JoinLayers(ss)

# sketch the cortical subset of the Visium HD dataset
DefaultAssay(ss) <- "Spatial"
ss <- FindVariableFeatures(ss)
ss <- SketchData(ss,ncells = 20000,method = "LeverageScore",sketched.assay = "sketch")

DefaultAssay(ss) <- "sketch"
ss <- ScaleData(ss)
ss <- RunPCA(ss, assay = "sketch", reduction.name = "pca.ss.sketch", verbose = T)
ss <- FindNeighbors(ss, reduction = "pca.ss.sketch", dims = 1:20)
ss <- RunUMAP(ss, reduction = "pca.ss.sketch", reduction.name = "umap.ss.sketch", return.model = T, dims = 1:20, verbose = T)

pdf(paste0("figurefolder/DimPlot_Lessthan2years.pdf"), height = 9, width = 13)
print(DimPlot(ss, reduction = "umap", group.by = c("ident", "orig.ident")))
dev.off()

pdf(paste0("figurefolder/SpatialDimPlot_Lessthan2years.pdf"), height = 9, width = 13)
print(SpatialDimPlot(ss))
dev.off()

pdf(paste0("figurefolder/SpatialFeaturePlot_", nfile,".pdf"), height = 9, width = 13)
print(SpatialFeaturePlot(ss, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")))
dev.off()

# load in the reference scRNA-seq dataset
load("H:/Backup_I_J_drive/ScottNess/ACC/integrated.RData")
integrated[["RNA"]] <- JoinLayers(integrated[["RNA"]])
Idents(integrated) <- "HPCAType"
counts <- integrated[["RNA"]]$counts
cluster <- as.factor(integrated$HPCAType)
nUMI <- integrated$nCount_RNA

#seurat_obj[[nfile]][["RNA"]] <- split(seurat_obj[[nfile]][["RNA"]], f = seurat_obj[[nfile]]$stim)

#Integration with scRNA-seq data (deconvolution)
#Robust Cell Type Decomposition, a computational approach to deconvolve spot-level data from spatial datasets, when provided with an scRNA-seq reference. RCTD has been shown to accurately annotate spatial data from a variety of technologies, including SLIDE-seq, Visium, and the 10x Xenium in-situ spatial platform. We observe good performance with Visium HD as well.

#To run RCTD, we first install the spacexr package from GitHub which implements RCTD
if (!requireNamespace("spacexr", quietly = TRUE)) {
  devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
}
library(spacexr)
# create the RCTD reference object
reference <- Reference(counts, cluster, nUMI)

counts_hd <- ss[["sketch"]]$counts
cortex_cells_hd <- colnames(ss[["sketch"]])
coords <- GetTissueCoordinates(ss)[cortex_cells_hd, 1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 28)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
# add results back to Seurat object
cortex <- AddMetaData(cortex, metadata = RCTD@results$results_df)

# project RCTD labels from sketched cortical cells to all cortical cells
cortex$first_type <- as.character(cortex$first_type)
cortex$first_type[is.na(cortex$first_type)] <- "Unknown"
cortex <- ProjectData(
  object = cortex,
  assay = "Spatial",
  full.reduction = "pca.cortex",
  sketched.assay = "sketch",
  sketched.reduction = "pca.cortex.sketch",
  umap.model = "umap.cortex.sketch",
  dims = 1:50,
  refdata = list(full_first_type = "first_type")
)

DefaultAssay(object) <- "Spatial"

# we only ran RCTD on the cortical cells
# set labels to all other cells as "Unknown"
object[[]][, "full_first_type"] <- "Unknown"
object$full_first_type[Cells(cortex)] <- cortex$full_first_type[Cells(cortex)]

Idents(object) <- "full_first_type"

# now we can spatially map the location of any scRNA-seq cell type
# start with Layered (starts with L), excitatory neurons in the cortex
cells <- CellsByIdentities(object)
excitatory_names <- sort(grep("^L.* CTX", names(cells), value = TRUE))
p <- SpatialDimPlot(object, cells.highlight = cells[excitatory_names], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T, ncol = 4)
p

``` 


for (mfile in files){
    
    #Highlight the spatial localization of a few clusters 
    #cells <- CellsByIdentities(ss, idents = levels(ss))
    cells <- CellsByIdentities(ss, idents = c(0,1,2,3,4))
    for (i in 1:5){
        pdf(paste0("figurefolder/Spatial_cluster", i,"_",nfile,".pdf"), height = 9, width = 13)
        print(SpatialDimPlot(ss,cells.highlight = cells[setdiff(names(cells)[i], "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend())
        dev.off()
    }
    
}
