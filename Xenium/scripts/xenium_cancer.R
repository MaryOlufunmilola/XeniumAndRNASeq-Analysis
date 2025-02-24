# -----------------------------------------------------------------------------------------------
# Load libraries
# -----------------------------------------------------------------------------------------------

PROJECTNAME <- ""
setwd(paste0("/home/ubuntu/", PROJECTNAME)

# Install and Load required libraries
# install.packages("pacman", repos='http://cran.us.r-project.org', quiet=TRUE)
# library(devtools)
# install_github('immunogenomics/presto')
# bio_pkgs <- c('SingleR','glmGamPoi','celldex') 
# BiocManager::install(bio_pkgs,force=TRUE)
# invisible(lapply(bio_pkgs, function(x) library(x, character.only=TRUE)))
      
library(pacman)
pacman::p_load(Polychrome, future, ggplot2, ggrepel, gridExtra, Seurat,scCustomize, UCell,presto,parallel,dplyr)             

xeniumData <- "data/"  # args[1:(length(args) - 1)]
resultsDir <- "results/" # args[length(args)]
if (!file.exists(resultsDir)){ dir.create(resultsDir) }  

options(Seurat.object.assay.version = "v5") 
no_cores <- availableCores() - 4 ##setting options for future
options(future.globals.maxSize = 230 * 1024^3, future.rng.onMisuse = 'ignore') 

n <- 400
set.seed(50)
mycolors = as.vector(createPalette(n,  c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', 
                                         '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', 
                                         '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'), target = "normal", M=50000))

# -----------------------------------------------------------------------------------------------
# Starting pre-processing
# -----------------------------------------------------------------------------------------------

#Use pattern that corresponds to data for analysis
DATAPATTERN <- ""
xfiles = list.files(xeniumData, pattern = DATAPATTERN)
xenium.obj = list()

# Loading Xenium dataset
for (j in xfiles) {
  seu <- ReadXenium(data.dir = paste0(xeniumData,j), outs = c("matrix", "microns"),type = c("centroids", "segmentations"))
  seu$segmentations$cell <- paste0(j, "_",seu$segmentations$cell)
  seu$centroids$cell <- paste0(j, "_",seu$centroids$cell)
  colnames(seu$matrix$`Gene Expression`) <- paste0(j, "_",colnames(seu$matrix$`Gene Expression`))
  colnames(seu$matrix$`Negative Control Probe`) <- paste0(j, "_",colnames(seu$matrix$`Negative Control Probe`))
  colnames(seu$matrix$`Negative Control Codeword`) <- paste0(j, "_",colnames(seu$matrix$`Negative Control Codeword`))
  xenium.obj[[j]] <- CreateSeuratObject(counts = seu$matrix[["Gene Expression"]], assay = "Xenium")
  xenium.obj[[j]][["ControlCodeword"]] <- CreateAssay5Object(counts = seu$matrix[["Negative Control Codeword"]])
  xenium.obj[[j]][["ControlProbe"]] <- CreateAssay5Object(counts = seu$matrix[["Negative Control Probe"]])
  xenium.obj[[j]][[paste0(j,"_fov")]] <- CreateFOV(
                                                  coords = list(centroids = CreateCentroids(seu$centroids),
                                                                segmentation = CreateSegmentation(seu$segmentations)),
                                                  type = c("segmentation", "centroids"),
                                                  molecules = seu$microns,
                                                  assay = "Xenium")
  xenium.obj[[j]]$sample <- j
 
# QC plot of genes per cell (nFeature_Xenium) and transcript counts per cell (nCount_Xenium).
  png(paste0(resultsDir,"VlnPlotbeforefilter_",j,".png"), height = 9, width = 13, units = 'in', res = 300)
  print(VlnPlot(xenium.obj[[j]], features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0,cols = mycolors,raster=FALSE))
  dev.off()

  # save sample as rds file
  saveRDS(xenium.obj[[j]], paste0(resultsDir,j,".rds"))

  # remove cells with less than 20 counts
  xenium.obj[[j]] <- subset(xenium.obj[[j]], subset = nCount_Xenium > 20)
  
  xenium.obj[[j]] <- RenameCells(xenium.obj[[j]], add.cell.id = j)
  
}

# -----------------------------------------------------------------------------------------------
# Downstream Analysis 1
# -----------------------------------------------------------------------------------------------
# Merge all 8 samples
xen_integrated <- merge(xenium.obj[[1]], y = xenium.obj[2:8])

# Save merged seurat object
saveRDS(xen_integrated, paste0(resultsDir,"xenium_merged.rds"))
rm(xenium.obj)

# Running SCTransform-based normalization
xen_integrated <- SCTransform(xen_integrated,assay = "Xenium")

# Running PCA,CLustering,UMAP reductions
xen_integrated <- RunPCA(xen_integrated, npcs = 30, features = rownames(xen_integrated)) 
xen_integrated <- FindNeighbors(xen_integrated, reduction = "pca", dims = 1:30)
xen_integrated <- FindClusters(xen_integrated, resolution = 0.3, cluster.name = "uninteg.clusters.r3")
xen_integrated <- RunUMAP(xen_integrated, dims = 1:30, reduction="pca", reduction.name = "umap.unintegrated")

# Integration of 8 samples
xen_integrated <- IntegrateLayers(object = xen_integrated, method = HarmonyIntegration, orig.reduction = "pca", 
                                    new.reduction = "harmony", assay = "SCT", verbose = TRUE)
# Clustering
xen_integrated <- FindNeighbors(xen_integrated, reduction = "harmony", dims = 1:30)
xen_integrated <- FindClusters(xen_integrated, resolution = 0.03, cluster.name = "integ.clusters.r3")

# Running UMAP reductions
xen_integrated <- RunUMAP(xen_integrated, dims = 1:30, reduction = "harmony", reduction.name = "umap")

# Save integrated object
saveRDS(xen_integrated, paste0(resultsDir,"xenium_integrated_final.rds"))

Idents(xen_integrated) <- xen_integrated$integ.clusters.r3
xen_integrated$ordered_clusters <- factor(xen_integrated$integ.clusters.r3, 
                                          levels = c("0","1","2","3","4","5","6","7","8","9","10",
                                                     "11","12","13","14","15","16","17","18","19",
                                                     "20","21","22","23","24","25"))
Idents(xen_integrated) <- xen_integrated$ordered_clusters

# Save the Elbow Plot as a PNG image with specified dimensions and resolution
png("ElbowPlot.png", height = 9, width = 13, units = 'in', res = 300)
print(ElbowPlot(xen_integrated, ndims=30,reduction="pca"))
dev.off() 

# Save the Dim Plots as a PNG image 
png("DimPlotSampleb4.png", height = 9, width = 13, units = 'in', res = 300)
print(DimPlot(xen_integrated, reduction= "umap.unintegrated", pt.size = 0.01,group.by = "sample", raster=FALSE, cols=mycolors) + ggtitle("Xenium Samples"))
dev.off() 

png("DimPlotClusterb4.png", height = 9, width = 13, units = 'in', res = 300)
print(DimPlot(xen_integrated, label = T, reduction= "umap.unintegrated", label.size= 1, pt.size = 0.01,group.by = "uninteg.clusters.r3", raster=FALSE, cols=mycolors) + ggtitle("Xenium Clusters"))
dev.off()

png("DimPlotSample.png", height = 9, width = 13, units = 'in', res = 300)
print(DimPlot(xen_integrated, reduction= "umap", pt.size = 0.01,group.by = "sample", raster=FALSE, cols=mycolors) + ggtitle("Xenium Clusters by Sample"))
dev.off()

png("DimPlotcluster.png", height = 9, width = 13, units = 'in', res = 300)
print(DimPlot(xen_integrated, label = T, reduction= "umap", label.size= 4, pt.size = 0.01, shuffle=T, 
                group.by = "ordered_clusters",raster=FALSE, cols=mycolors) + ggtitle("Xenium Clusters"))
dev.off() 

png("DimPlotcluster_bySample.png", height = 9, width = 20, units = 'in', res = 300)
print(DimPlot(xen_integrated, label = T, reduction= "umap", label.size= 3, ncol=4,pt.size = 0.01,shuffle=T, group.by = "ordered_clusters", 
                split.by = "sample", raster=FALSE, cols=mycolors) + ggtitle("Xenium Clusters by Sample"))
dev.off() 

#presto::wilcoxauc compares cells in each cluster to all other cells in the dataset. 
xen_integrated = PrepSCTFindMarkers(xen_integrated, assay = "SCT", verbose = TRUE)

all_markers <- presto::wilcoxauc(xen_integrated,  seurat_assay = "SCT", assay = "data", group_by = "ordered_clusters") %>%
                    Add_Pct_Diff(pct.1_name = "pct_in",  pct.2_name = "pct_out")
write.csv(all_markers, paste0(resultsDir,"xenium_integrated_markergenes.csv"), quote = F)

# Find markers expressed in greater than 50% of target population
top10 <- presto::top_markers(all_markers, n = 10, auc_min = 0.5, pval_max = 0.05, padj_max = 0.05,pct_in_min = 20)
top_markers_auc <- unique(c(na.omit(unlist(top10[,-1],use.names=F))))

png("top7markerspercluster_dotplotauc.png", width = 20, height = 9, units = 'in', res = 300)
print(DotPlot_scCustom(xen_integrated, features = top_markers_auc, group.by = "integ.clusters.r3",x_lab_rotate = TRUE) + 
              theme(axis.text.x=element_text(size=4),axis.text.y=element_text(size=6)))
dev.off()

png("Heatmap_topmarkersauc.png", width = 20, height = 9, units = 'in', res = 300)
print(DoHeatmap(xen_integrated, features = top_markers_auc, group.by = "integ.clusters.r3",
                cells = Random_Cells_Downsample(xen_integrated, num_cells = 100, allow_lower = T))
     + theme(axis.text.y=element_text(size=3)))
dev.off()

# Create Violin Plots
top_features <- c("DCN", "FLT1", "CXCL8", "IDO1", "CXCL13","EPCAM","EGR3", "MMP12", "KIT", "NT5E","CCL21", "IGHM", "IGKC", "MKI67", "RGS5",
                  "S100A9",  "IGHG2", "IL7R", "CD8A", "CD4","CD3E", "CD19", "MS4A1", "CD68", "VSIG4","KLRD1", "CD247", "CLEC10A", "FCER1A")
png("VlnPlots.png", height = 9, width = 13, units = 'in', res = 300)
print(Stacked_VlnPlot(xen_integrated, features = top_features, x_lab_rotate = TRUE,group.by = "ordered_clusters"))
dev.off() 


