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
DimPlot(choroidCN.combined, reduction = "umap", split.by = "batch")

table(choroidCN.combined$batch)

DefaultAssay(choroidCN.combined) <- "RNA"

FeaturePlot(choroidCN.combined, features = "nCount_RNA", min.cutoff = "q9")

#Quantitate cell numMembers
data <- choroidCN.combined
Idents(data) <- "type"
controls <- subset(data, ident = "Normal")
paps <- subset(data, ident = "CPP")
cars <- subset(data, ident = "CPC")
total <- table(data$seurat_clusters)
Control <- table(controls$seurat_clusters)
CPP <- table(paps$seurat_clusters)
CPC <- table(cars$seurat_clusters)
cellcounts <- cbind(Control, CPP, CPC, total)
write.table(cellcounts, file = "ChoroidTumorCellCounts.csv", sep = ",")

#marker expression objects
marker_survey <- c("OTX2", "AQP1", "PECAM1", "FLT1", "CD163", "CD74", "LAMA2", "FN1", "COL1A2", "CENPE", "ASPM", "NRXN1", "SYT1", "GRIK1", "GRIN2A", "CHGB", "CCK",  "CD247", "CD96")
insitu <- c("CD163", "CD274", "CD96", "FLT1", "ANGPT2", "ANGPT1", "FN1", "CSPG4", "TOP2A", "MKI67" )
WNTgenes1 <- c("WNT1", "WNT2", "WNT2B", "WNT3", "WNT3A", "WNT4", "WNT5A", "WNT5B", "WNT6", "WNT7A", "WNT7B", "WNT8A")
WNTgenes2 <- c("WNT8B", "WNT10A", "WNT10B", "WNT11", "WNT14", "WNT15", "WNT16")
proliferation <- c("CENPF", "ASPM", "CENPE", "TOP2A", "MKI67")
CILIAgenes <- c("FOXJ1", "DNAAF1", "KIF24", "MAK", "DRC1", "TCTE1", "IQCD", "DNAH2", "DNALI1", "DNAI1", "TEKT", "TEKT2", "TEKT4")
MacActivation1 <- c("ARG1", "MGL2", "TMEM26", "RNASE2A", "MRC1", "EGR2", "FLT1", "CHIL3", "CLEC10A", "ATP6V0D2")
MacActivation2 <- c("FAM198b", "MATK", "SOCS2", "ITGB3", "OCSTAMP", "PTGS1", "S100A4", "CLEC7A", "PLXDC2", "HBEGF", "CCL24", "EMP1", "PDCD1LG2", "OLFM1", "UBE2C")
MacActivationRnD1 <- c("CD80", "CD86", "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL11", "CCL19", "CCL20", "CD36", "CD68", "PTGS2")
MacActivationRnD2 <- c("CX3CL1", "CXCL1", "CXCL5", "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CXCL16", "FCGR3A", "FCGR2A", "IFNG", "IFNGR1")
MacActivationRnD3 <- c("IL1B", "IL6", "CXCL8", "IL12A", "IL15", "IL17A", "IL18", "IL23A", "IRF5")
TcellMemory <- c("CD44", "CD62L", "CD127", "CXCR3", "PDCD1", "CD274", "PDCD1LG2")
MHCs <- c("HLA-DRB5", "HLA-E", "HLA-DRB1", "HLA-DRA", "HLA-DQB1", "HLA-DPB1", "HLA-DPA1", "HLA-DMB", "HLA-C")
MPvsMG <- c("CD14", "CD163", "CD16", "CCR2", "APOE", "P2RY12", "TREM2", "GPR34", "CXCR3", "ITGAM", "ITGAX")
CSCmarkers <- c("ALDH1A1", "CD44", "CD47", "CEACAM1", "BMI1", "CD24", "PROM1", "GJA1", "THY1", "CXCR4")
hippoTargets1 <- c("ATAD2", "BMP4", "CTGF", "CYR61", "FGF2", "FRZD1", "IRS1", "PTF2", "SHISA9", "TILL1")
YAP1c7 <- c("LINC00662", "ACTG1", "RPLP0", "RPL6", "RPL4", "RPS19", "RPS3", "TOP2A", "CTGF", "RPS15", "RPL7", "DIAPH3", "RPL41", "RPL10", "TUBB", "RPS11", "ID1", "RPS6", "RPL13A", "RPL37A", "RPL15", "RPS27A", "RPL28", "ACAT1", "IGF1R", "ZEB1", "SGK1")
Endothelial_progenitors <- c("CDH5", "CD34", "PROM1", "CXCR4", "ETV2", "MCAM","VCAM1", "NT5E", "PTPRC", "KIT", "TEK", "FLT1", "KDR")
c15 <- c("KCNIP4", "PCDH9", "ADGRB3", "LRRTM4", "NCAM2", "ROBO2", "CADM2", "GRIK2", "DSCAM", "ADARB2", "DPP10", "GFAP", "GRIA4", "NRXN1", "NCAM1", "GRM5", "KCNH7", "SYT1", "MAP2")


Tumor_clusters <- c(0,1,)
#Find Conserved markers in clusters
markers <- vector()
UC <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20 )
for (clusters in UC) {
  marker <- FindConservedMarkers(choroidCN.combined, ident.1 = clusters,  grouping.var = "status", verbose = TRUE)
  file_name <- paste0(clusters,"CNallmarkers800.csv")
  write.table(marker, file = file_name, sep = ",")
  print(file_name)
  markers <- c(markers, marker)
}
save(markers,file='choroidCN_combinedMarkers-all.Robj')

#Find distinguishing markers in similar cell types
MarcrophageSubtypes <- FindMarkers(choroidCN.combined, ident.1 = 8,  ident.2 = 6, verbose = TRUE)
bloodV1 <- FindMarkers(choroidCN.combined, ident.1 = 5, ident.2 = 9, verbose = FALSE)
bloodV2 <- FindMarkers(choroidCN.combined, ident.1 = 5, ident.2 = 12, verbose = FALSE)
bloodV3 <- FindMarkers(choroidCN.combined, ident.1 = 9, ident.2 = 12, verbose = FALSE)
##Carcinoma enriched clusters
for (cluster in c(2,3,10,14)) {
  cancergenes <- FindMarkers(choroidCN.combined, ident.1 = cluster, ident.2 = 0, verbose = FALSE)
  file_name <- paste0(cluster,"Carcinomagenes.csv")
  write.table(cancergenes, file = file_name, sep = ",")
}

##Papilloma enriched clusters
for (cluster in c(7,11,16,18)) {
  cancergenes <- FindMarkers(choroidCN.combined, ident.1 = cluster, ident.2 = 0, verbose = FALSE)
  file_name <- paste0(cluster,"Papillomagenes.csv")
  write.table(cancergenes, file = file_name, sep = ",")
}

##Other tumor specific clusters
for (cluster in c(13,17,20)) {
  cancergenes <- FindMarkers(choroidCN.combined, ident.1 = cluster, ident.2 = 0, verbose = FALSE)
  file_name <- paste0(cluster,"Tumorgenes.csv")
  write.table(cancergenes, file = file_name, sep = ",")
}

data <- choroidCN.combined
Idents(data) <- "seurat_clusters"
Carclusters23 <- subset(data, ident = c("0", "1", "2", "3"))
VlnPlot(Carclusters23, features = c("PROM1", "CDH12", "JUN", "RPS27A", "SFRP1", "TGFBR2", "WNT2B", "WNT5B", "WNT9A" ), split.by = "batch", group.by = "seurat_clusters", pt.size = 0, combine = FALSE)

MacClusters <- subset(data, iden = c("6", "8"))
VlnPlot(MacClusters, features = MacActivation1, split.by = "batch", group.by = "seurat_clusters", pt.size = 0, combine = FALSE)

VlnPlot(choroidCN.combined, features = c("WNT5B","OTX2", "CD163", "CD247", "CD79A", "MKI67", "VWF",  "PECAM1"), split.by = "type", group.by = "seurat_clusters", pt.size = 0, combine = FALSE)


# Annotate cell types
choroidCN.renamed <- choroidCN.combined
choroidCN.renamed <- RenameIdents(choroidCN.combined, `0` = "Epithelial", `1` = "Epithelial", `2` = "Epithelial", 
                                   `3` = "Epithelial", `4` = "Epithelial", `5` = "Pericyte1", `6` = "Macrophage1", `7` = "PapEnrich", `8` = "Macrophage2", `9` = "Endothelial", 
                                   `10` = "Proliferative", `11` = "PapEnrich", '12' = "Pericyte2", '13' = "CarEnrich", '14' = "CarEnrich", '15' = "NeuroEndocrine", '16' = "PapEnrich", '17' = "CarEnrich", '18' = "PapEnrich", '19' = "T cells", '20' = "CarEnrich")

DimPlot(choroidCN.renamed, reduction = "umap", label = TRUE)


#confirm cell type assignment
choroidCN.renamed$celltype <- Idents(choroidCN.renamed)
VlnPlot(choroidCN.renamed, features = marker_survey, split.by = "type", group.by = "celltype", pt.size = 0, combine = FALSE)
FeaturePlot(object = choroidCN.renamed, features = marker_survey, cols.use = c("grey", "blue"), 
            reduction.use = "umap")
#compare cell type specific gene expression
## get tumor type DEGs
###Epithelial Cells

CPT.combined <- choroidCN.renamed
CPT.combined$celltype.type <- paste(Idents(CPT.combined), CPT.combined$type, sep = "_")
Idents(CPT.combined) <- "celltype.type"
Epi.response <- FindMarkers(CPT.combined, ident.1 = "Epithelial_CPC", ident.2 = "Epithelial_CPP", verbose = FALSE)
head(Epi.response, n = 15)
write.table(Epi.response, file = "EpithelialDEGS.csv", sep = ",")

Epi.responseCC <- FindMarkers(CPT.combined, ident.1 = "Epithelial_CPC", ident.2 = "Epithelial_Normal", verbose = FALSE)
head(Epi.responseCC, n = 15)
write.table(Epi.responseCC, file = "EpithelialCCDEGS.csv", sep = ",")

Epi.responsePC <- FindMarkers(CPT.combined, ident.1 = "Epithelial_CPP", ident.2 = "Epithelial_Normal", verbose = FALSE)
head(Epi.responseCC, n = 15)
write.table(Epi.responsePC, file = "EpithelialPCDEGS.csv", sep = ",")

CarNorSigDegs <- filter(Epi.responseCC, abs(avg_logFC) > 1)
celltype <- "Epithelial"
comparison <- "Car_v_Norm"
DEG_file <- paste0(celltype, comparison, "DEGS")
  
genes <- as.numeric(CarNorSigDegs$p_val_adj)
names(genes) <- rownames(CarNorSigDegs)
GOdata <- new("topGOdata",
              ontology = "BP", # use biological process ontology
              allGenes = genes,
              geneSelectionFun = function(x)(x < .05),
              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")

resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

fischTable <- GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
write.table(fischTable, file = paste0(DEG_file, "fishT.csv"), sep = ",")


g <- ggplot(fischTable, aes(x=reorder(Term, Annotated), y=Annotated , label=Expected)) + 
  geom_bar(stat='identity', aes(fill=Term), width=.5,position="dodge")  + theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = "none") + 
  labs(title= paste0("Top 20 enriched GO terms ", celltype, " ", comparison), x="GO term", y="Number Genes (Fisher's statistic)") + coord_flip() + geom_text(aes(label=paste0("(", Fisher, ")")), position=position_dodge(width=0.1), vjust=0.25, hjust=0)
g

  
  
CPT.combined <- choroidCN.renamed
CPT.combined$celltype.status <- paste(Idents(CPT.combined), CPT.combined$status, sep = "_")
Idents(CPT.combined) <- "celltype.status"
Epi.responseTumor <- FindMarkers(CPT.combined, ident.1 = "Epithelial_Tumor", ident.2 = "Epithelial_Normal", verbose = FALSE)
head(Epi.responseTumor, n = 15)
write.table(Epi.responseTumor, file = "EpithelialTumorDEGS.csv", sep = ",")
####Tumor specific epithelial clusters
TumorSpec <- subset(data, ident = c("7", "10", "11", "13", "14", "16", "17", "20","0")) #0=control
VlnPlot(TumorSpec, features = c("WNT5B","PROM1"), split.by = "batch", group.by = "seurat_clusters", pt.size = 0, combine = FALSE)



###Macrophages
#Cell number ANOVA
group_by(MacTable, Type) %>%
  summarise(
    count = n(),
    mean = mean(percentM2, na.rm = TRUE),
    sd = sd(percentM2, na.rm = TRUE)
  )

library("ggpubr")
ggboxplot(MacTable, x = "Type", y = "percentM2", 
          color = "Type", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          order = c("Normal", "Pap", "Car"),
          ylab = "%M2 macrophage", xlab = "Tissue Type")

# Compute the analysis of variance
res.aov <- aov(percentM2 ~ Type, data = MacTable)
# Summary of the analysis
summary(res.aov)

TukeyHSD(res.aov)


CPT.combined <- choroidCN.combined
CPT.combined$celltype.type <- paste(Idents(CPT.combined), CPT.combined$type, sep = "_")
Idents(CPT.combined) <- "celltype.type"
Mac.response6CN <- FindMarkers(CPT.combined, ident.1 = "6_CPC", ident.2 = "6_Normal", verbose = FALSE)
head(Mac.response6CN, n = 15)
write.table(Mac.response6CN, file = "Macrophage6CNDEGS.csv", sep = ",")

Mac.response6CP <- FindMarkers(CPT.combined, ident.1 = "6_CPC", ident.2 = "6_CPP", verbose = FALSE)
head(Mac.response6CP, n = 15)
write.table(Mac.response6CP, file = "Macrophage6CPDEGS.csv", sep = ",")

Mac.response6PN <- FindMarkers(CPT.combined, ident.1 = "6_CPP", ident.2 = "6_Normal", verbose = FALSE)
head(Mac.response6PN, n = 15)
write.table(Mac.responsePC, file = "Macrophage6PNDEGS.csv", sep = ",")

Mac.response8CN <- FindMarkers(CPT.combined, ident.1 = "8_CPC", ident.2 = "8_Normal", verbose = FALSE)
head(Mac.response6CN, n = 15)
write.table(Mac.response8CN, file = "Macrophage8CNDEGS.csv", sep = ",")

Mac.response8CP <- FindMarkers(CPT.combined, ident.1 = "8_CPC", ident.2 = "8_CPP", verbose = FALSE)
head(Mac.response8CP, n = 15)
write.table(Mac.response8CP, file = "Macrophage8CPDEGS.csv", sep = ",")

Mac.response8PN <- FindMarkers(CPT.combined, ident.1 = "8_CPP", ident.2 = "8_Normal", verbose = FALSE)
head(Mac.response8PN, n = 15)
write.table(Mac.response8PN, file = "Macrophage8PNDEGS.csv", sep = ",")


T.responseCP <- FindMarkers(CPT.combined, ident.1 = "19_CPC", ident.2 = "19_CPP", verbose = FALSE)
head(T.responseCP, n = 15)
write.table(T.responseCP, file = "TcellDEGSCP.csv", sep = ",")

T.responseCN <- FindMarkers(CPT.combined, ident.1 = "19_CPC", ident.2 = "19_Normal", verbose = FALSE)
head(T.responseCN, n = 15)
write.table(T.responseCN, file = "TcellDEGSCN.csv", sep = ",")

T.responsePN <- FindMarkers(CPT.combined, ident.1 = "19_CPP", ident.2 = "19_Normal", verbose = FALSE)
head(T.responsePN, n = 15)
write.table(T.responsePN, file = "TcellDEGSPN.csv", sep = ",")


CPT.combined <- choroidCN.renamed
CPT.combined$celltype.status <- paste(Idents(CPT.combined), CPT.combined$status, sep = "_")
Idents(CPT.combined) <- "celltype.status"
MAC.responseTumor <- FindMarkers(CPT.combined, ident.1 = "Macrophage_Tumor", ident.2 = "Macrophage_Normal", verbose = FALSE)
head(MAC.responseTumor, n = 15)
write.table(MAC.responseTumor, file = "MacrophageTumorDEGS.csv", sep = ",")

# Go analysis

## Create topGOdata object
Cluster_no <- 7
cancer_type <- "Carcinoma"
DEG_file <- paste0("C", Cluster_no, cancer_type, "N.response.csv")

genes <- as.numeric(C7PN.response$p_val_adj)
names(genes) <- rownames(C7PN.response)
GOdata <- new("topGOdata",
              ontology = "BP", # use biological process ontology
              allGenes = genes,
              geneSelectionFun = function(x)(x < .05),
              annot = annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol")

resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

fischTable <- GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
write.table(fischTable, file = paste0(DEG_file, "fishT.csv"), sep = ",")


g <- ggplot(fischTable, aes(x=reorder(Term, Annotated), y=Annotated , label=Expected)) + 
  geom_bar(stat='identity', aes(fill=Term), width=.5,position="dodge")  + theme(plot.title = element_text(hjust = 0.5, size = 14), legend.position = "none") + 
  labs(title= paste0("Top 20 enriched GO terms cluster ", Cluster_no, " ", cancer_type), x="GO term", y="Number Genes (Fisher's statistic)") + coord_flip() + geom_text(aes(label=paste0("(", Fisher, ")")), position=position_dodge(width=0.1), vjust=0.25, hjust=0)
g
