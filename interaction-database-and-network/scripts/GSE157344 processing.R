library(Seurat)
library(stringr)
library(dplyr)
library(Cairo)
library(Matrix)
sample <- list.files(pattern = ".h5$")
sceList = lapply(sample, function(x){
  print(x)
  p <- strsplit(x,split = "_")[[1]][1]
  a=Read10X_h5(x,use.names = TRUE)
  b = try(a$`Gene Expression`,silent = TRUE)
  if("try-error" %in% class(b)){
    sce <- CreateSeuratObject(a,project = x, min.features=200)
  }else{sce <- CreateSeuratObject(a$`Gene Expression`,project = p, min.features=200)}
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
sce <- FindNeighbors(sce, dims = 1:30)
sce <- FindClusters(sce, resolution = 0.2)
DimPlot(sce, reduction = "umap", group.by = "seurat_clusters", raster = FALSE)

inn.markers <- c("EPCAM", # epithelium
                 "PTPRC","CD3G","CD3E","CD79A" # immune
)
DotPlot(sce, features = inn.markers, group.by = "integrated_snn_res.0.2", assay = "RNA") + RotatedAxis()
sce$cell.type <- ifelse(sce$integrated_snn_res.0.2 %in% 5, "Airway epithelial","immune")


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
immu <- FindClusters(immu, resolution = 1.5)
DotPlot(immu, features = immu.markers, group.by = "integrated_snn_res.1.5", assay = "RNA") + RotatedAxis()
immu$cell.type[immu$integrated_snn_res.1 %in% c(25,26)] <- "B cell"
immu$cell.type[immu$integrated_snn_res.1 %in% c(5,23)] <- "Plasma cell"
immu$cell.type[immu$integrated_snn_res.1 %in% c(11,12,14,22)] <- "T cell"
immu$cell.type[immu$integrated_snn_res.1 %in% 20] <- "NK cell"
immu$cell.type[immu$integrated_snn_res.1 %in% c(0,1,2,4,7,10,15,21,24)] <- "Neutrophil"
immu$cell.type[immu$integrated_snn_res.1 %in% c(6,13,17,19)] <- "Macrophage"
immu$cell.type[immu$integrated_snn_res.1 %in% 3] <- "Dendritic cell"
immu$cell.type[immu$integrated_snn_res.1 %in% c(8,9,16,18)] <- "Monocyte"
DimPlot(immu, reduction = "tsne", group.by = "cell.type")
sce$cell.type[sce$cell.type == "immune"] <- immu$cell.type

sp <- as.character(unique(sce$orig.ident))
for(i in sp){
  sample <- subset(sce, orig.ident == i)
  saveRDS(sample,paste(i,".rds",sep = ""))
}

allsample <- list.files(pattern =  ".rds$")[-22:-23]
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
