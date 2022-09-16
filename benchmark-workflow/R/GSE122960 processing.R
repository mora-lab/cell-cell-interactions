rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(stringr)
library(clustree)
fs = list.files(pattern = '.h5')
fs

# Create Seurat data list
sceList = lapply(fs, function(x){
  print(x)
  a=Read10X_h5(x)
  p=str_split(x,'_',simplify = T)
  p=paste(p[2],p[3],sep = "_")
  sce <- CreateSeuratObject( a ,project = p )
  sce
})
for (i in 1:length(sceList)) {
  sceList[[i]] <- NormalizeData(sceList[[i]], verbose = FALSE)
  sceList[[i]] <- FindVariableFeatures(sceList[[i]], selection.method = "vst")
}

# Integrate data
sceList
sce.anchors <- FindIntegrationAnchors(object.list = sceList, dims = 1:30)
sce.integrated <- IntegrateData(anchorset = sce.anchors, dims = 1:30)

# Normalization, clustering
sce.integrated <- ScaleData(sce.integrated, verbose = FALSE)
sce.integrated <- RunPCA(sce.integrated, npcs = 30, verbose = FALSE)
ElbowPlot(sce.integrated, ndims = 30)
sce.integrated <- RunUMAP(sce.integrated, reduction = "pca", dims = 1:30)
sce.integrated <- RunTSNE(sce.integrated, reduction = "pca", dims = 1:30)
DimPlot(sce.integrated, reduction = "tsne", group.by = "orig.ident")

meta <- data.frame(matrix(data=NA, nrow=ncol(sce.integrated),ncol=2))
rownames(meta) <- colnames(sce.integrated)
colnames(meta) <- c("ID", "Class")
meta[,1] <- sce.integrated$orig.ident
meta[grep("Donor",meta[,1]),"Class"] <- "Donor"
meta[grep("IPF",meta[,1]),"Class"] <- "IPF"
sce.integrated$Class <- meta[,"Class"]

sce.integrated <- FindNeighbors(sce.integrated, dims = 1:30)
sce.integrated <- FindClusters(sce.integrated, resolution = c(seq(0,1.2,0.2)))
clustree(sce.integrated@meta.data, prefix = "integrated_snn_res.")
sce.integrated <- FindClusters(sce.integrated, resolution = 0.8)
DimPlot(sce.integrated, reduction = "tsne", group.by = "seurat_clusters")

# Annotation
all.markers <- c("IGKC","IGLC2","IGHG3","IGHG4","IGLC3", # Plasma
                 "APOC1","C1QB","FABP4","C1QA","APOE", # Macro
                 "SFTPC","SFTPA1","SFTPA2","PGC","NAPSA", # AT2
                 "SCGB1A1","SCGB3A1","SCGB3A2","BPIFB1","MMP7", # Club
                 "MS4A1","LTB","CD79A","CD79B","AL928768.3", # B
                 "CCL5","IL32","GNLY","NKG7","CCL4", # T/NK
                 "S100A8","S100A9","FCN1","S100A12","VCAN", # Mono
                 "CCL17","S100B","HLA-DPB1","RGS1","MS4A6A", # Dendritic
                 "S100A2","KRT15","KRT17","MMP1","KRT5", # Basal
                 "TPSB2","TPSAB1","CPA3","MS4A2","KIT", # Mast
                 "SPARCL1","CCL21","CLDN5","VWF","IGFBP7", # Endothelial/Lymphatic
                 "A2M","DCN","MGP","COL1A2","COL3A1","ACTA2", # Fibro
                 "AGER","CAV1","EMP2","NTM","CEACAM6", # AT1
                 "CAPS","C9orf24","TPPP3","C20orf85","TMEM190" # Ciliated
)
levels(Idents(sce.integrated))
new.cluster.ids <- c("Macrophage","Macrophage","AT2","AT2","AT2",
                     "Plasma","Club/Basal","Dendritic cell","Monocyte","Macrophage",
                     "Ciliated","Macrophage","AT1","Tcell","Endothelial",
                     "Macrophage","Unknown","B cell","Macrophage","Fibroblast",
                     "Mast","AT2","B cell","Endothelial","AT2","Macrophage","AT2")
names(new.cluster.ids) <- levels(sce.integrated)
sce.integrated <- RenameIdents(sce.integrated, new.cluster.ids)
table(Idents(sce.integrated))
sce.integrated$cell.type <- Idents(sce.integrated)
DimPlot(sce.integrated, reduction = "umap", group.by = "cell.type")

# Save data for each IPF sample
IPFsample <- subset(sce.integrated, Class == "IPF")
cell <- colnames(IPFsample[,Idents(IPFsample) %in% c("Plasma","Club/Basal","Ciliated","Unknown","B cell", "Dendritic cell")])
IPFsample <- subset(IPFsample, cells = cell, invert= T) 
sampledata <- unique(IPFsample$orig.ident)
for(x in sampledata){
  sub <- subset(IPFsample, orig.ident == x)
  saveRDS(sub,paste("GSE122960/", x, ".rds", sep = ""))
}
