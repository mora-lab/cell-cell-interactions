library(Seurat)
library(stringr)
library(dplyr)
library(Cairo)
library(Matrix)
sample <- list.files()[]
sceList = lapply(sample, function(x){
  print(x)
  a=Read10X(x)
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
sce <- FindClusters(sce, resolution = 0.2)

inn.markers <- c("EPCAM", # epithelium
                 "PECAM1", # Endothelial
                 "COL1A1","PDGFRA", # Fibroblast
                 "CNN1","ACTA2","TAGLN","RGS5","DES", # Smooth Muscle
                 "PTPRC","CD3G","CD3E","CD79A" # immune
)
DotPlot(sce, features = inn.markers, group.by = "integrated_snn_res.0.2", assay = "RNA") + RotatedAxis()
sm <- "9"
epi <- c("5","12","13","17")
endo <- c("6","10","14","26")
fibro <- "3"
sce$cell.type <- ifelse(sce$integrated_snn_res.0.2 %in% sm, "Smooth muscle",
                        ifelse(sce$integrated_snn_res.0.2 %in% epi, "epi", 
                               ifelse(sce$integrated_snn_res.0.2 %in% endo, "Endothelial",
                                      ifelse(sce$integrated_snn_res.0.2 %in% fibro, "Fibroblast","immune"))))

epi.markers <- c("AGER","PDPN","CLIC5","CLDN18","FSTL3", # Alveolar Epithelial Type 1 Cell
                 "SFTPB","SFTPA1","SFTPD","MUC1","ETV5" # Alveolar Epithelial Type 2 Cell
)
epi <- subset(sce,cell.type == "epi")
epi <- ScaleData(epi, verbose = FALSE)
epi <- RunPCA(epi, npcs = 30, verbose = FALSE)
epi <- RunUMAP(epi, reduction = "pca", dims = 1:30)
epi <- RunTSNE(epi, reduction = "pca", dims = 1:30)
epi <- FindNeighbors(epi, dims = 1:30)
epi <- FindClusters(epi, resolution = 0.2)
DotPlot(epi, features = epi.markers, group.by = "integrated_snn_res.0.2", assay = "RNA") + RotatedAxis()
at1 <- "1"
at2 <- c("0","2","4","5","6","7","9","10","11")
epi$cell.type <- ifelse(epi$integrated_snn_res.0.2 %in% at1, "AT1",
                        ifelse(epi$integrated_snn_res.0.2 %in% at2, "AT2","Airway epithelial"))
DimPlot(epi, reduction = "tsne", group.by = "cell.type")
sce$cell.type[sce$cell.type == "epi"] <- epi$cell.type


immu.markers <- c("CD79A","CD24","MS4A1","CD19", # B Cell
                  "CD27", "SLAMF7", # Plasma Cell
                  "CD3E","CD3D","CD3G","TRBC2", # T Cell
                  "KLRD1","NKG7","TRDC","NCR1", # Natural Killer Cell
                  "S100A8","S100A9","LCN2", "MMP9", "CEACAM8", "FCGR3B", "CXCR2", # Neutrophil
                  "MS4A2","CPA3","TPSAB1", # Mast Cell
                  "MARCO","MSR1","MRC1","CD68", # Macrophage
                  "CD83","CD86","ZBTB46", # Dendritic Cell
                  "CD14","APOBEC3A","CFP","VCAN","FCN1" # Monocyte
)
immu <- subset(sce,cell.type == "immune")
immu <- ScaleData(immu, verbose = FALSE)
immu <- RunPCA(immu, npcs = 30, verbose = FALSE)
immu <- RunUMAP(immu, reduction = "pca", dims = 1:30)
immu <- RunTSNE(immu, reduction = "pca", dims = 1:30)
immu <- FindNeighbors(immu, dims = 1:30)
immu <- FindClusters(immu, resolution = 0.5)
DotPlot(immu, features = immu.markers, group.by = "integrated_snn_res.0.5", assay = "RNA") + RotatedAxis()
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(12,15)] <- "B cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(19,20)] <- "Plasma cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(1,2,3,16,17)] <- "T cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(5,8)] <- "NK cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(11,18)] <- "Neutrophil"
immu$cell.type[immu$integrated_snn_res.0.5 %in% 6] <- "Mast cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(4,10)] <- "Macrophage"
immu$cell.type[immu$integrated_snn_res.0.5 %in% 13] <- "Dendritic cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(0,9)] <- "Monocyte"
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(7,14)] <- "Unknown"
DimPlot(immu, reduction = "umap", group.by = "cell.type")
sce$cell.type[sce$cell.type == "immune"] <- immu$cell.type
sce <- subset(sce,cell.type == "Unknown",invert=TRUE)

sp <- as.character(unique(sce$orig.ident))
for(i in sp){
  sample <- subset(sce, orig.ident == i)
  saveRDS(sample,paste(i,".rds",sep = ""))
}

allsample <- list.files(pattern =  ".rds$")[-4]
for (sample in allsample){
  seu <- readRDS(sample)
  CellPhoneDB_process(seu,"cpdb/")
}

CellPhoneDB_process = function(seu,dir){
  # save data
  sample <- unique(seu$orig.ident)
  dir <- paste(dir, sample, sep = "")
  dir.create(dir)
  writeMM(seu@assays$RNA@counts, file = paste(dir,"/matrix.mtx",sep=""))
  write(x = rownames(seu@assays$RNA@counts), file = paste(dir,"/features.tsv",sep=""))
  write(x = colnames(seu@assays$RNA@counts), file = paste(dir,"/barcodes.tsv",sep=""))
  
  # save meta
  seu@meta.data$Cell = rownames(seu@meta.data)
  df <- seu@meta.data[, c('Cell', 'cell.type')]
  write.table(df, file = paste(dir,"/meta.tsv",sep=""), sep = '\t', quote = F, row.names = F)
}
