library(Seurat)
library(stringr)
library(clustree)
library(dplyr)

sample <- list.files("raw")
data_dir <- paste("raw", sample, sep="/")

# Create Seurat data
sceList = lapply(data_dir, function(x){
  print(x)
  a=Read10X(x)
  p=str_split(x,"/",simplify = T)
  p=p[2]
  sce <- CreateSeuratObject( a ,project = p, min.features=200)
  sce
})

for (i in 1:length(sceList)) {
  sceList[[i]][["percent.mt"]] <- PercentageFeatureSet(sceList[[i]], pattern = "^MT-")
  sceList[[i]] <- subset(sceList[[i]], subset = percent.mt < 35)
  sceList[[i]] <- NormalizeData(sceList[[i]], verbose = FALSE)
  sceList[[i]] <- FindVariableFeatures(sceList[[i]], selection.method = "vst")
}
sceList

# Integrate data
sce.anchors <- FindIntegrationAnchors(object.list = sceList, dims = 1:30, reduction = "cca")
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
meta[grep("NOR",meta[,1]),"Class"] <- "Donor"
meta[grep("IPF",meta[,1]),"Class"] <- "IPF"
sce.integrated$Class <- meta[,"Class"]

sce.integrated <- FindNeighbors(sce.integrated, dims = 1:30)
sce.integrated <- FindClusters(sce.integrated, resolution = 0.8)
DimPlot(sce.integrated, reduction = "tsne", group.by = "seurat_clusters")

# Annotation
all.markers <- c("SCGB1A1","SCGB3A2","AGER", "SFTPC", "MUC5B", "FOXJ1", "KRT5",
                 "COL1A1","COL1A2","PDGFRA","DES","ACTG2","VWF","RGS5","LYVE1",
                 "AIF1","CD163","CD8A","CD79A","MS4A1","IGKC","IGHA1","IGHG3","MZB1",
                 "CD1C","GNLY","NKG7","GZMB","PRF1","CST7","CD3D","TRAC","CCL3","CCL4","CCL5","TPSAB1","CPA3","MS4A2","FCN1"
)
levels(Idents(sce.integrated))
new.cluster.ids <- c("Macrophage","Macrophage","Macrophage","Unknown","Unknown",
                     "Tcell","Fibroblast","Endothelial","Monocyte","AT2",
                     "Ciliated cell","Dendritic cell","Club/Goblet cell","Pericyte","AT1",
                     "Unknown","NK cell","Mast","Unknown","Fibroblast",
                     "Smooth muscle cell","B cell","Ciliated cell","Basal cell","Dendritic cell",
                     "Unknown","Lymphatic endothelial","Monocyte","Unknown","AT1")
names(new.cluster.ids) <- levels(sce.integrated)
sce.integrated <- RenameIdents(sce.integrated, new.cluster.ids)
table(Idents(sce.integrated))
sce.integrated$cell.type <- Idents(sce.integrated)
DimPlot(sce.integrated, reduction = "tsne", group.by = "cell.type")

# # Save data for each IPF sample
IPFsample <- subset(sce.integrated, Class == "IPF")
cell <- colnames(IPFsample[,Idents(IPFsample) %in% c("Lymphatic endothelial","Pericyte","NK cell","Ciliated cell","Dendritic cell",
                                               "Club/Goblet cell","Basal cell","Unknown","B cell","Smooth muscle cell")])
IPFsample <- subset(IPFsample, cells = cell, invert= T)
sampledata <- unique(IPFsample$orig.ident)
for(x in sampledata){
  sub <- subset(IPFsample, ori.ident == x)
  saveRDS(sub,paste("GSE128033/", x, ".rds", sep = ""))
}