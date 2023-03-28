library(Seurat)
library(stringr)
library(dplyr)
library(Cairo)
library(Matrix)
options(bitmapType='cairo')
mtx <- fread("GSE171541_UMI_count.txt", header = T)
mtx <- as.data.frame(mtx)
mtx <- mtx[!duplicated(mtx[,1]),]
rownames(mtx) <- mtx[,1]
mtx[,1] <- NULL
mtx <- as(mtx,"dgCMatrix")
sce <- CreateSeuratObject(counts = mtx)
# remove healthy smoker
sce <- subset(sce,orig.ident == "ctl.o1" | orig.ident == "ctl.o2", invert=TRUE)
saveRDS(sce,"sce_raw.rds")
sce <- NormalizeData(sce, verbose = FALSE)
sce <- FindVariableFeatures(sce, selection.method = "vst")
sce <- ScaleData(sce, verbose = FALSE)
sce <- RunPCA(sce, npcs = 30, verbose = FALSE)
DimPlot(sce, reduction = "pca", group.by = "orig.ident")
ElbowPlot(sce, ndims = 30)
sce<- RunUMAP(sce, reduction = "pca", dims = 1:30)
sce <- RunTSNE(sce, reduction = "pca", dims = 1:30)
sce <- FindNeighbors(sce, dims = 1:30)
sce <- FindClusters(sce, resolution = 0.2)
DimPlot(sce, reduction = "umap", group.by = "seurat_clusters")

inn.markers <- c("EPCAM", # epithelium
                 "PECAM1", # Endothelial
                 "COL1A1","PDGFRA", # Fibroblast
                 "CNN1","ACTA2","TAGLN","RGS5","DES", # Smooth Muscle
                 "PTPRC","CD3G","CD3E","CD79A" # immune
)
DotPlot(sce, features = inn.markers, group.by = "RNA_snn_res.0.2") + RotatedAxis()
imm <- c("1","3","4","7","9","13","14","17")
epi <- c("2","5","10","11")
endo <- c("0","6","8","18")
fibro <- c("12","15")
sce$cell.type <- ifelse(sce$RNA_snn_res.0.2 %in% imm, "immune",
                        ifelse(sce$RNA_snn_res.0.2 %in% epi, "epi", 
                               ifelse(sce$RNA_snn_res.0.2 %in% endo, "Endothelial",
                                      ifelse(sce$RNA_snn_res.0.2 %in% fibro, "Fibroblast","Smooth muscle"))))
DimPlot(sce, reduction = "tsne", group.by = "cell.type")

epi.markers <- c("AGER","PDPN","CLIC5","CLDN18","FSTL3", # Alveolar Epithelial Type 1 Cell
                 "SFTPB","SFTPA1","SFTPD","MUC1","ETV5" # Alveolar Epithelial Type 2 Cell
)
epi <- subset(sce,cell.type == "epi")
epi <- NormalizeData(epi, verbose = FALSE)
epi <- FindVariableFeatures(epi, selection.method = "vst")
epi <- ScaleData(epi, verbose = FALSE)
epi <- RunPCA(epi, npcs = 30, verbose = FALSE)
epi <- RunUMAP(epi, reduction = "pca", dims = 1:30)
epi <- RunTSNE(epi, reduction = "pca", dims = 1:30)
epi <- FindNeighbors(epi, dims = 1:30)
epi <- FindClusters(epi, resolution = 0.2)
DotPlot(epi, features = epi.markers, group.by = "RNA_snn_res.0.2") + RotatedAxis()
at1 <- c("4")
at2 <- c("0","2","3","6","8","9")
epi$cell.type <- ifelse(epi$RNA_snn_res.0.2 %in% at1, "AT1",
                        ifelse(epi$RNA_snn_res.0.2 %in% at2, "AT2", "Airway epithelial"))
DimPlot(epi, reduction = "umap", group.by = "cell.type")
sce$cell.type[sce$cell.type == "epi"] <- epi$cell.type


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
immu <- NormalizeData(immu, verbose = FALSE)
immu <- FindVariableFeatures(immu, selection.method = "vst")
immu <- ScaleData(immu, verbose = FALSE)
immu <- RunPCA(immu, npcs = 30, verbose = FALSE)
immu <- RunUMAP(immu, reduction = "pca", dims = 1:30)
immu <- RunTSNE(immu, reduction = "pca", dims = 1:30)
immu <- FindNeighbors(immu, dims = 1:30)
immu <- FindClusters(immu, resolution = 0.5)
DotPlot(immu, features = immu.markers, group.by = "RNA_snn_res.0.5") + RotatedAxis()
immu$cell.type[immu$RNA_snn_res.0.5 %in% c(11,13)] <- "B cell"
immu$cell.type[immu$RNA_snn_res.0.5 %in% c(14,16)] <- "Plasma cell"
immu$cell.type[immu$RNA_snn_res.0.5 %in% c(1,2,10)] <- "T cell"
immu$cell.type[immu$RNA_snn_res.0.5 %in% c(0,6,15)] <- "NK cell"
immu$cell.type[immu$RNA_snn_res.0.5 %in% 9] <- "Neutrophil"
immu$cell.type[immu$RNA_snn_res.0.5 %in% c(4,19)] <- "Mast cell"
immu$cell.type[immu$RNA_snn_res.0.5 %in% c(5,12,18)] <- "Macrophage"
immu$cell.type[immu$RNA_snn_res.0.5 %in% c(3,17)] <- "Dendritic cell"
immu$cell.type[immu$RNA_snn_res.0.5 %in% c(7,8)] <- "Monocyte"
DimPlot(immu, reduction = "tsne", group.by = "cell.type")
sce$cell.type[sce$cell.type == "immune"] <- immu$cell.type

sp <- as.character(unique(sce$orig.ident))
for(i in sp){
  sample <- subset(sce, orig.ident == i)
  saveRDS(sample,paste(i,".rds",sep = ""))
}

allsample <- c(paste("healthy/",list.files("healthy/"),sep=""), 
               paste("copd/",list.files("copd/"),sep=""))
dir.create("healthy/cpdb")
for (sample in allsample){
  seu <- readRDS(sample)
  CellPhoneDB_process(seu,"healthy/cpdb/")
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

