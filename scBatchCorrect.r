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
tumors <- c("C1", "C3", "T1", "T2", "T3", "T4", "T5", "T6", "T7", "T8", "T9", "T10", "T11", "T12")
tumor_types <- vector(mode = "list", length = 14)
tumor_types <- c("Normal", "Normal", "CPP", "CPP", "CPP", "aCPP", "CPC", "CPC", "CPC", "CPC", "CPC", "CPP", "CPP", "CPC")
names(tumor_types) <- tumors
status_list <- vector(mode = "list", length = 14)
status_list <- c('Normal', 'Normal', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor', 'Tumor','Tumor', 'Tumor')
names(status_list) <- tumors
#tumors <- c("C1", "C3", "T2", "T11")
tumors <- c("C1", "C3", "T2", "T3", "T4", "T6", "T7", "T8", "T12")
data.folder = "/icgc/dkfzlsdf/analysis/OE0519_projects/chptumor/filteredmatrices/"
for (tumor in tumors){
  file_loc <- paste0(data.folder,tumor,"/filtered_feature_bc_matrix")
  data <- Read10X(data.dir = file_loc)
  # Initialize the Seurat object with the raw (non-normalized data).
  object <- CreateSeuratObject(counts = data, project = "BatchCorrect", min.cells = 3, min.features = 200)
  #filter cells
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 10)
  object
  object$batch <- tumor
  object$type <- tumor_types[tumor]
  object$status <- status_list[tumor]
  samples <- c(samples, object)
  print(tumor)
}




#filter genes

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
  sample.celNorm <- FindVariableFeatures(sample.celNorm, selection.method = "vst", nfeatures = 2000)
  samples.cNorm <- c(samples.cNorm, sample.celNorm)
  S.no <- S.no+1
  print(S.no)
  print(length(samples.cNorm))
}

nsamples <- vector()
for (sample in samples){
  sample <- subset(sample, subset = nFeature_RNA > 400)
  sample <- NormalizeData(sample, verbose = TRUE)
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 2000)
  nsamples <- c(nsamples, sample)
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
# file_names <- tumors
# for (names in file_names){
#   fname <- paste0(names,".Robj")
#   load(file=fname)
#   print(fname)
# }

file_no <- 1
for (object in samples.cNorm){
  fname <- fname <- paste0(file_names[file_no],"CN.Robj")
  save(object,file=fname)
  print(fname)
  file_no <- file_no+1
}

#batch correction
choroid.anchors <- FindIntegrationAnchors(object.list = nsamples, dims = 1:20)
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
choroidCN.combined <- IntegrateData(anchorset = choroid.anchors, dims = 1:20)

choroidCN.combined <- ScaleData(choroidCN.combined, verbose = FALSE)
choroidCN.combined <- RunPCA(choroidCN.combined, npcs = 30, verbose = TRUE)

DefaultAssay(choroidCN.combined) <- "integrated"
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

#Find Conserved markers in clusters
markers <- vector()
for (clusters in levels(choroidCN.combined$seurat_clusters)) {
  marker <- FindConservedMarkers(choroidCN.combined, ident.1 = clusters, grouping.var = "type", verbose = TRUE)
  file_name <- paste0(clusters,"markers.csv")
  write.csv(marker, file = file_name, row.names=TRUE)
  print(file_name)
  markers <- c(markers, marker)
}
save(markers,file='choroidCN_combinedMarkers.Robj')


#Epithelial
FeaturePlot(choroidCN.combined, features = c("KCNQ1", "TTR", "OTX2", "TJP1", "AQP1", "EGFL7", "KIF18B", "PARD3", "MUC1", 
                                             "PRKCI"), min.cutoff = "q9")

#Endothelial
FeaturePlot(choroidCN.combined, features = c("NOTCH4", "EGFL7", "VWF", "DSP", "EFNB2", "BTNL9", "DLL4", "COL4", "FLT1", "FLRT2", "PECAM1", "VWF", "COL2A1", 
                                             "COL3A2", "CDH2", "SNAI1", "EPHA7", "SLC24A3", "SULT1E1", "SLC4A4", "PDZRN4", "THSD4", "FN1", "TWIST1", "COL4A1", "COL4A2"), min.cutoff = "q9")
#Mesenchymal
FeaturePlot(choroidCN.combined, features = c("CARMN", "SLC6A1", "EBF2", "CDH6", "DACH1", "COL1A2"), min.cutoff = "q9")

#Proliferative?
FeaturePlot(choroidCN.combined, features = c("ASPM", "KIF18B", "KIF14", "TOP2A", "KIF2C", "TACC3", "NCAPG", "C2orf48", "CENPF", "MKI67"), min.cutoff = "q9")

#TCell
FeaturePlot(choroidCN.combined, features = c("CD247", "SLFN12L", "CD2", "LCK", "CD96", "SCML4", "ITK", "THEMIS", "NGK"), min.cutoff = "q9")
#Bcell
FeaturePlot(choroidCN.combined, features = c("CD79A", "BTLAL", "FCRL2", "LCK", "IGHD", "FCRL5", "CD22", "MS4A1", "FCRL1"), min.cutoff = "q9")
#Macrophages&Monocytes
FeaturePlot(choroidCN.combined, features = c("CD163", "SLC11A1", "PIK3R5",  "RGS1", "DOCK2", "LAPTM5", "C3", "CSF3R" ,"TRPM2", "SLCO2B1"), min.cutoff = "q9")

choroidCN.combined <- RenameIdents(choroidCN.combined, `0` = "Epithelial", `1` = "Epithelial", `2` = "Epithelial", 
                                   `3` = "Epithelial", `4` = "Immune", `5` = "Epithelial", `6` = "Epithelial", `7` = "Epithelial", `8` = "Epithelial", `9` = "Epithelial", 
                                   `10` = "Epithelial", `11` = "Epithelial", `12` = "Endothelial", `13` = "Mesenchymal", `14` = "Mesenchymal", `15` = "Immune" )

#compare cell type specific gene expression
## get tumor type DEGs
CPT.combined <- choroidCN.combined
CPT.combined$celltype <- Idents(CPT.combined)
Idents(CPT.combined) <- "celltype.type"
Epi.response <- FindMarkers(CPT.combined, ident.1 = "Epithelial_CPC", ident.2 = "Epithelial_CPP", verbose = FALSE)
head(Epi.response, n = 15)
write.table(Epi.response, file = "EpithelialDEGS.csv", sep = ",")


e.cells <- subset(choroidCN.combined, idents = "Epithelial")
Idents(e.cells) <- "type"
avg.e.cells <- log1p(AverageExpression(e.cells, verbose = FALSE)$RNA)
avg.e.cells$gene <- rownames(avg.e.cells)

genes.to.label = c("WNT5B", "HIF1A", "YAP1", "MKI67", "WNT2B", "VEGF", "IL6", "CXCL10", "CCL8")
p1 <- ggplot(avg.e.cells, aes(CPP, CPC)) + geom_point() + ggtitle("Epithelial Cells")

#Cell type specific marker expression
proliferative <- c("MKI67", "TOP2A", "ASPM", "TACC3")
epithelial <- c("AQP1", "OTX2", "PARD3", "KCNQ1", "TTR")
VlnPlot(object = CPT.combined, features = proliferative, split.by = "type", group.by = "celltype")
VlnPlot(object = CPT.combined, features = epithelial, split.by = "type", group.by = "celltype")

library(data.table)
data_to_write_out <- as.data.frame(as.matrix(choroidCN.combined[[RNA]]))
data_to_write_out <- as.data.frame(as.matrix(choroidCN.combined@assays$RNA))
fwrite(x = data_to_write_out, row.names = TRUE, file = "outfile.txt")
