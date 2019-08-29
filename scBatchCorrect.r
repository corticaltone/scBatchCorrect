# scRNAseq analysis of human normal and tumor ChP 
library(scater) # BioConductor
library(loomR) # BioConductor
library(dplyr)
library(Seurat)
library(DropletUtils) # BioConductor
library(scran)
library(ggplot2)
library(cowplot)

# Load the ChP datasets
samples <- vector()
tumors <- c("C1", "C3", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T11", "T12")
tumor_types <- vector(mode = "list", length = 12)
tumor_types <- c("Normal", "Normal", "CPP", "CPP", "CPP", "aCPP", "CPC", "CPC", "CPC", "CPC", "CPP", "CPC")
names(tumor_types) <- tumors
status_list <- vector(mode = "list", length = 12)
status_list <- c('Normal', 'Normal', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor')
names(status_list) <- tumors
tumors <- c("T9", "T10")
data.folder = "/icgc/dkfzlsdf/analysis/OE0519_projects/chptumor/filteredmatrices/"
for (tumor in tumors){
  file_loc <- paste0(data.folder,tumor,"/filtered_feature_bc_matrix")
  data <- Read10X(data.dir = file_loc)
  # Initialize the Seurat object with the raw (non-normalized data).
  object <- CreateSeuratObject(counts = data, project = "BatchCorrect", min.cells = 3, min.features = 200)
  object
  object$batch <- tumor
  object$type <- tumor_types[tumor]
  samples <- c(samples, object)
  object$status <- status_list[tumor]
  print(tumor)
}

#filter
for (sample in samples){
  sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")
  sample <- subset(sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
}

#make sce
samples.sce <- vector()
for (sample in samples){
  sce <- as.SingleCellExperiment(sample)
  samples.sce <- c(samples.sce, sce)
}

#normalize
samples.cNorm <- vector()
S.no <- 0
for (sample in samples.sce){
  print(sample)
  clusters <- quickCluster(sample, min.size=100)
  sample <- computeSumFactors(sample, cluster=clusters)
  sample <- normalize(sample)
  sample.celNorm <- as.Seurat(sample, counts = "counts", data = "logcounts")
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
  samples.cNorm <- c(samples.cNorm, sample.celNorm)
  S.no <- S.no+1
  print(S.no)
  print(length(samples.cNorm))
}

for (sample in samples){
  sample <- subset(sample, subset = nFeature_RNA > 400)
  sample <- NormalizeData(sample, verbose = TRUE)
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
}



#Save Sueurat objects

file_names <- tumors
file_no <- 1
for (object in samples){
  fname <- paste0(file_names[file_no],".Robj")
  save(object,file=fname)
  print(fname)
  file_no <- file_no+1
}

# Reentry point 
file_names <- tumors
file_no <- 1
for (object in samples){
  fname <- paste0(file_names[file_no],".Robj")
  load(object,file=fname)
  print(fname)
  file_no <- file_no+1
}

file_no <- 1
for (object in samples.cNorm){
  fname <- fname <- paste0(file_names[file_no],"CN.Robj")
  save(object,file=fname)
  print(fname)
  file_no <- file_no+1
}

#batch correction
choroid.anchors <- FindIntegrationAnchors(object.list = samples, dims = 1:20)
choroid.combined <- IntegrateData(anchorset = choroid.anchors, dims = 1:20)
DefaultAssay(choroid.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
choroid.combined <- ScaleData(choroid.combined, verbose = FALSE)
choroid.combined <- RunPCA(choroid.combined, npcs = 30, verbose = TRUE)

# t-SNE and Clustering
choroid.combined <- RunUMAP(choroid.combined, reduction = "pca", dims = 1:20)
choroid.combined <- FindNeighbors(choroid.combined, reduction = "pca", dims = 1:20)
choroid.combined <- FindClusters(choroid.combined, resolution = 0.5)
save(choroid.combined,file='choroid_combined.Robj')

#Repeat for cell normalized Data
choroidCN.anchors <- FindIntegrationAnchors(object.list = samples.cNorm, dims = 1:20)
choroidCN.combined <- IntegrateData(anchorset = choroidCN.anchors, dims = 1:20)
DefaultAssay(choroidCN.combined) <- "integrated"
choroidCN.combined <- ScaleData(choroidCN.combined, verbose = FALSE)
choroidCN.combined <- RunPCA(choroidCN.combined, npcs = 30, verbose = FALSE)

choroidCN.combined <- RunUMAP(choroidCN.combined, reduction = "pca", dims = 1:20)
choroidCN.combined <- FindNeighbors(choroidCN.combined, reduction = "pca", dims = 1:20)
choroidCN.combined <- FindClusters(choroidCN.combined, resolution = 0.5)

#Save clustered objects
save(choroid.combined,file='choroid_combined.Robj')
save(choroidCN.combined,file='choroidCN_combined.Robj')

# Visualization
p1 <- DimPlot(choroidCN.combined, reduction = "umap", group.by = "batch")
p2 <- DimPlot(choroidCN.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(choroidCN.combined, reduction = "umap", split.by = "batch")

DefaultAssay(choroidCN.combined) <- "RNA"
FeaturePlot(choroidCN.combined, features = c("CD163", "TTR", "OTX2", "AQP1", "SULF1", "PECAM1", "VWF", "COL2A1", 
                                           "COL3A2", "PPBP"), min.cutoff = "q9")

FeaturePlot(choroidCN.combined, features = c("CD163", "TTR", "CDH1", "DSP", "TJP1", "OTX2", "AQP1", "SULF1", "FLT1", "CDH5", "PECAM1", "VWF", "COL2A1", 
                                             "COL3A2", "CDH2", "SNAI1", "EPHA7", "SLC24A3", "SULT1E1", "SLC4A4", "PDZRN4", "THSD4", "FN1", "TWIST1"), min.cutoff = "q9")

choroid4.markers <- FindConservedMarkers(choroidCN.combined, ident.1 = 4, grouping.var = "batch", verbose = FALSE)
head(choroid4.markers)

choroid6.markers <- FindConservedMarkers(choroidCN.combined, ident.1 = 6, grouping.var = "batch", verbose = FALSE)
head(choroid6.markers)

choroidCN.combined <- RenameIdents(choroidCN.combined, `0` = "Epithelial", `1` = "Epithelial", `2` = "Epithelial", 
                                `3` = "Epithelial", `4` = "Immune", `5` = "Epithelial", `6` = "Epithelial", `7` = "Epithelial", `8` = "Epithelial", `9` = "Epithelial", 
                                `10` = "Epithelial", `11` = "Epithelial", `12` = "Endothelial", `13` = "Mesenchymal", `14` = "Mesenchymal", `15` = "Immune" )

DimPlot(choroidCN.combined, label = TRUE)

i.cells <- subset(choroidCN.combined, idents = "Immune")
Idents(i.cells) <- "batch"
avg.i.cells <- log1p(AverageExpression(i.cells, verbose = FALSE)$RNA)
avg.i.cells$gene <- rownames(avg.i.cells)

genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p1 <- ggplot(avg.i.cells, aes(C1, T1)) + geom_point() + ggtitle("CD4 Naive T Cells")

T3 <- subset(choroidCN.combined, subset = batch == "T3")

library(data.table)
data_to_write_out <- as.data.frame(as.matrix(choroidCN.combined[[RNA]]))
data_to_write_out <- as.data.frame(as.matrix(choroidCN.combined@assays$RNA))
fwrite(x = data_to_write_out, row.names = TRUE, file = "outfile.txt")
