library(Seurat)
library(stringr)
library(dplyr)
library(Cairo)
library(Matrix)
sample <- list.files()[-14]
sceList = lapply(sample, function(x){
  print(x)
  a=Read10X_h5(x)
  sce <- CreateSeuratObject(a,project = x, min.features=200)
  sce
})

for (i in 1:length(sceList)) {
  sceList[[i]][["percent.mt"]] <- PercentageFeatureSet(sceList[[i]], pattern = "^MT-")
  sceList[[i]] <- subset(sceList[[i]], subset = percent.mt < 35)
  sceList[[i]] <- NormalizeData(sceList[[i]], verbose = FALSE)
  sceList[[i]] <- FindVariableFeatures(sceList[[i]], selection.method = "vst")
}
sceList

sce.anchors <- FindIntegrationAnchors(object.list = sceList, dims = 1:30, reduction = "cca")
sce <- IntegrateData(anchorset = sce.anchors, dims = 1:30)
saveRDS(sce,"sce_raw.rds")
sce <- ScaleData(sce, verbose = FALSE)
sce <- RunPCA(sce, npcs = 30, verbose = FALSE)
DimPlot(sce, reduction = "pca", group.by = "orig.ident")
ElbowPlot(sce, ndims = 30)
sce<- RunUMAP(sce, reduction = "pca", dims = 1:30)
sce <- RunTSNE(sce, reduction = "pca", dims = 1:30)
sce <- FindNeighbors(sce, dims = 1:30)
sce <- FindClusters(sce, resolution = 0.8)

inn.markers <- c("EPCAM", # epithelium
                 "PECAM1", # Endothelial
                 "COL1A1","PDGFRA", # Fibroblast
                 "CNN1","ACTA2","TAGLN","RGS5","DES", # Smooth Muscle
                 "PTPRC","CD3G","CD3E","CD79A" # immune
)
DotPlot(sce, features = inn.markers, group.by = "integrated_snn_res.0.8", assay = "RNA") + RotatedAxis()
sm <- 14
epi <- c(4,7,12,22,26,33)
endo <- c(9,19,24,29,31)
fibro <- 15
sce$cell.type <- ifelse(sce$integrated_snn_res.0.8 %in% sm, "Smooth muscle",
                        ifelse(sce$integrated_snn_res.0.8 %in% epi, "epi", 
                               ifelse(sce$integrated_snn_res.0.8 %in% endo, "Endothelial",
                                      ifelse(sce$integrated_snn_res.0.8 %in% fibro, "Fibroblast","immune"))))
DimPlot(sce, reduction = "tsne", group.by = "cell.type")

epi.markers <- c("AGER","PDPN","CLIC5","CLDN18","FSTL3", # Alveolar Epithelial Type 1 Cell
                 "SFTPB","SFTPA1","SFTPD","MUC1","ETV5", # Alveolar Epithelial Type 2 Cell
                 "TTR","CD44","FABP5","TFF3","IGFBP5","EGFR"
)
epi <- subset(sce,cell.type == "epi")
epih <- subset(epi,orig.ident == "p028t.h5"|orig.ident == "p039t.h5"|orig.ident == "p044t.h5",invert=TRUE)
epit<-  subset(epi,orig.ident == "p028t.h5"|orig.ident == "p039t.h5"|orig.ident == "p044t.h5")
epih <- ScaleData(epih, verbose = FALSE)
epih <- RunPCA(epih, npcs = 30, verbose = FALSE)
epih <- RunUMAP(epih, reduction = "pca", dims = 1:30)
epih <- RunTSNE(epih, reduction = "pca", dims = 1:30)
epih <- FindNeighbors(epih, dims = 1:30)
epih <- FindClusters(epih, resolution = 0.2)
DotPlot(epih, features = epi.markers, group.by = "integrated_snn_res.0.2", assay = "RNA") + RotatedAxis()
at1 <- 2
at2 <- 0
epih$cell.type <- ifelse(epih$integrated_snn_res.0.2 %in% at1, "AT1",
                        ifelse(epih$integrated_snn_res.0.2 %in% at2, "AT2", "Airway epithelial"))
epit <- ScaleData(epit, verbose = FALSE)
epit <- RunPCA(epit, npcs = 30, verbose = FALSE)
epit <- RunUMAP(epit, reduction = "pca", dims = 1:30)
epit <- RunTSNE(epit, reduction = "pca", dims = 1:30)
epit <- FindNeighbors(epit, dims = 1:30)
epit <- FindClusters(epit, resolution = 0.2)
DotPlot(epit, features = epi.markers, group.by = "integrated_snn_res.0.2", assay = "RNA") + RotatedAxis()
epit$cell.type <- "Tumor cell"

sce$cell.type[sce$cell.type == "epi" & !(sce$orig.ident == "p028t.h5"|sce$orig.ident == "p039t.h5"|sce$orig.ident == "p044t.h5")] <- epih$cell.type
sce$cell.type[sce$cell.type == "epi" & (sce$orig.ident == "p028t.h5"|sce$orig.ident == "p039t.h5"|sce$orig.ident == "p044t.h5")] <- epit$cell.type


immu.markers <- c("CD79A","CD24","MS4A1","CD19", # B Cell
                  "CD27", "SLAMF7", # Plasma Cell
                  "CD3E","CD3D","CD3G","TRBC2", # T Cell
                  "KLRD1","NKG7","TRDC","NCR1", # Natural Killer Cell
                  "S100A8","S100A9","FPR1","CSF3R","PTGS2", # Neutrophil
                  "MS4A2","CPA3","TPSAB1", # Mast Cell
                  "MARCO","MSR1","MRC1","CD68", # Macrophage
                  "CD83","CD86","ZBTB46", # Dendritic Cell
                  "CD14","APOBEC3A","CFP","MMP19","VCAN","FCN1" # Monocyte
)
immu <- subset(sce,cell.type == "immune")
immu <- ScaleData(immu, verbose = FALSE)
immu <- RunPCA(immu, npcs = 30, verbose = FALSE)
immu <- RunUMAP(immu, reduction = "pca", dims = 1:30)
immu <- RunTSNE(immu, reduction = "pca", dims = 1:30)
immu <- FindNeighbors(immu, dims = 1:30)
immu <- FindClusters(immu, resolution = 0.5)
DotPlot(immu, features = immu.markers, group.by = "integrated_snn_res.0.5", assay = "RNA") + RotatedAxis()
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(13,16)] <- "B cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(19,21,22)] <- "Plasma cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(1,4,5,15,17)] <- "T cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% 3] <- "NK cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% 6] <- "Neutrophil"
immu$cell.type[immu$integrated_snn_res.0.5 %in% 11] <- "Mast cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(0,7,10,14,18,20)] <- "Macrophage"
immu$cell.type[immu$integrated_snn_res.0.5 %in% 8] <- "Dendritic cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(2,9,12)] <- "Monocyte"
DimPlot(immu, reduction = "umap", group.by = "cell.type")
sce$cell.type[sce$cell.type == "immune"] <- immu$cell.type

sp <- as.character(unique(sce$orig.ident))
for(i in sp){
  sample <- subset(sce, orig.ident == i)
  n <- strsplit(i,split = "\\.")[[1]][1]
  saveRDS(sample,paste(n,".rds",sep = ""))
}

allsample <- list.files(pattern = ".rds$")[1:13]
dir.create("cpdb")
for (sample in allsample){
  seu <- readRDS(sample)
  CellPhoneDB_process(seu,"cpdb/")
}

CellPhoneDB_process = function(seu,dir){
  # save data
  sample <- unique(seu$orig.ident)
  n <- strsplit(sample,split = "\\.")[[1]][1]
  dir <- paste(dir, n, sep = "")
  dir.create(dir)
  writeMM(seu@assays$RNA@counts, file = paste(dir,"/matrix.mtx",sep=""))
  write(x = rownames(seu@assays$RNA@counts), file = paste(dir,"/features.tsv",sep=""))
  write(x = colnames(seu@assays$RNA@counts), file = paste(dir,"/barcodes.tsv",sep=""))
  
  # save meta
  seu@meta.data$Cell = rownames(seu@meta.data)
  df <- seu@meta.data[, c('Cell', 'cell.type')]
  write.table(df, file = paste(dir,"/meta.tsv",sep=""), sep = '\t', quote = F, row.names = F)
}
