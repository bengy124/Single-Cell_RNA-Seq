## ----Load libraries-----------------------------------------------------------
library(rmarkdown)
library(Seurat)
library(tidyverse)
library(dplyr)
library(hdf5r)

## ----Loading data-------------------------------------------------------------
# Loading Hodgkin's Lymphoma tumor dataset
braintumor.data<- Read10X_h5(filename = 'D:/My Drive/Projects/R/Single-Cell RNA-Seq/data/Brain_Tumor_3p_raw_feature_bc_matrix.h5')

# Initialize the Seurat object with the raw (non-normalized data)
braintumor <- CreateSeuratObject(counts = braintumor.data, project = "Hodgkin's Lymphoma", min.cells = 3, min.features = 200)
braintumor


## ----QC-----------------------------------------------------------------------
#Mitochondrial QC metrics
braintumor[["percent.mt"]] <- PercentageFeatureSet(braintumor, pattern = "^MT-")


## ----Plots--------------------------------------------------------------------
# Violin plot
VlnPlot(braintumor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
view(braintumor@meta.data)

# Scatter plot to visualize feature-feature relationships
plot1 <- FeatureScatter(braintumor, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(braintumor, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2


## ----Filtering----------------------------------------------------------------
braintumor <- subset(braintumor, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20)


## ----Normalization------------------------------------------------------------
# Normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result.
braintumor <- NormalizeData(braintumor)


## ----Feature selection--------------------------------------------------------
# Identifying variable features
braintumor <- FindVariableFeatures(braintumor, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(braintumor), 10)

# plot variable features with and without labels
plot3 <- VariableFeaturePlot(braintumor)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = T)
plot3
plot4

## ----Scaling------------------------------------------------------------------
all.genes <- rownames(braintumor)
braintumor <- ScaleData(braintumor, features = all.genes)

## ----PCA----------------------------------------------------------------------
braintumor <- RunPCA(braintumor, features = VariableFeatures(object = braintumor))

# Checking if cells cluster in PCA
DimPlot(braintumor, reduction = "pca") + NoLegend()

## ----Heatmap------------------------------------------------------------------
# Both cells and features are ordered according to their PCA scores
DimHeatmap(braintumor, dims = 1:3, cells = 500, balanced = T)

## ----Elbow plot---------------------------------------------------------------
# Elbow plot: a ranking of principle components based on the percentage of variance explained by each one 
ElbowPlot(braintumor)

## ----Clustering---------------------------------------------------------------
braintumor <- FindNeighbors(braintumor, dims = 1:15)
braintumor <- FindClusters(braintumor, resolution = c(0.1,0.3, 0.5, 0.7, 1))

# Visualization
DimPlot(braintumor, group.by = "RNA_snn_res.0.1", label = T)

# Setting identity of clusters
Idents(braintumor) <- "RNA_snn_res.0.1"


## ----Non-linear dimensionality reduciton--------------------------------------
braintumor <- RunUMAP(braintumor, dims = 1:15)
DimPlot(braintumor, reduction = "umap")

## ----Specific clusters--------------------------------------------------------
# Finding all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(braintumor, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# Finding markers for every cluster compared to all remaining cells, reporting only the significant ones
braintumor.markers <- FindAllMarkers(braintumor, only.pos = T)
braintumor.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Plotting raw counts
VlnPlot(braintumor, features =c("NOVA1", "CD14", "LYZ", "NKG7") , slot = "counts", log = T)

FeaturePlot(braintumor, features = c("CD14", "NOVA1", "FCGR3A", "LYZ", "CD8A", "IL7R"))

## ----Cell typing--------------------------------------------------------------
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono")
names(new.cluster.ids) <- levels(braintumor)
braintumor <- RenameIdents(braintumor, new.cluster.ids)
DimPlot(braintumor, reduction = "umap", label = T, pt.size = 0.5) + NoLegend()

