library(Seurat)
library(stringr)
library(dplyr)
library(Cairo)
library(Matrix)
library(data.table)
options(bitmapType='cairo')
cells <- read.table("head.txt")
group <- read.csv("GSE136246_human_cell_groupings.csv")
tumor <- c()
for(i in 1:41649){
  g <- strsplit(colnames(group)[i+1],split = "\\.")
  if(g[[1]][1] %in% c("C2","C10","C13")){
    tumor <- c(tumor, cells[1,i])
  }
}
sample <- list.files(pattern = ".tsv.gz$")
sceList = lapply(sample, function(x){
  print(x)
  m <- strsplit(x,split = "_")[[1]][2]
  n <- strsplit(m,split = "\\.")[[1]]
  if(length(n) == 5){
    p <- paste(n[1],n[2],sep = ".")
  }else{p <- n[1]}
  a=fread(x, header = TRUE)
  a=t(a)
  barcode <- a[1,]
  colnames(a) <- barcode
  a <- a[-1,]
  a <- as.data.frame(a)
  sce <- CreateSeuratObject(a,project = p, min.features=200)
  sce
})
saveRDS(sceList,"scelist.rds")
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
Idents(sce) <- newcol
sce$cell.type <- "NULL"
for(i in 1:65989){
  if(Idents(sce)[i] %in% tumor){
    sce$cell.type[i] <- "Tumor cell"}
}
subsce <- subset(sce,cell.type == "Tumor cell", invert=TRUE)


subsce <- ScaleData(subsce, verbose = FALSE)
subsce <- RunPCA(subsce, npcs = 30, verbose = FALSE)
subsce <- FindNeighbors(subsce, dims = 1:30)
subsce <- FindClusters(subsce, resolution = 1.0)
inn.markers <- c("EPCAM", # epithelium
                 "PECAM1", # Endothelial
                 "COL1A1","PDGFRA", # Fibroblast
                 "CNN1","ACTA2","TAGLN","RGS5","DES", # Smooth Muscle
                 "PTPRC","CD3G","CD3E","CD79A" # immune
)
DotPlot(subsce, features = inn.markers, group.by = "integrated_snn_res.1",assay = "RNA") + RotatedAxis()
fibrosm<- 19
epi <- c(6,7,10,32,33)
endo <- c(2,4,9,22,25)
subsce$cell.type <- ifelse(subsce$integrated_snn_res.1 %in% epi, "epi", 
                               ifelse(subsce$integrated_snn_res.1 %in% endo, "Endothelial",
                                      ifelse(subsce$integrated_snn_res.1 %in% fibrosm, "Fibroblast and Sm","imm")))


fism <- subset(subsce,cell.type=="Fibroblast and Sm")
fism <- ScaleData(fism, verbose = FALSE)
fism <- RunPCA(fism, npcs = 30, verbose = FALSE)
fism <- RunUMAP(fism, reduction = "pca", dims = 1:30)
fism <- RunTSNE(fism, reduction = "pca", dims = 1:30)
fism <- FindNeighbors(fism, dims = 1:30)
fism <- FindClusters(fism, resolution = 0.2)
DotPlot(fism, features = inn.markers, group.by = "integrated_snn_res.0.2", assay= "RNA") + RotatedAxis()
DotPlot(fism, features = immu.markers, group.by = "integrated_snn_res.0.2", assay= "RNA") + RotatedAxis()
fism$cell.type <- ifelse(fism$integrated_snn_res.0.2 %in% 0,"Fibroblast",
                         ifelse(fism$integrated_snn_res.0.2 %in% c(1,2),"Smooth muscle","Macrophage"))
subsce$cell.type[subsce$cell.type == "Fibroblast and Sm"]<- fism$cell.type


epi.markers <- c("AGER","PDPN","CLIC5","CLDN18","FSTL3", # Alveolar Epithelial Type 1 Cell
                 "SFTPB","SFTPA1","SFTPD","MUC1","ETV5" # Alveolar Epithelial Type 2 Cell
)
epi <- subset(subsce,cell.type == "epi")
epi <- ScaleData(epi, verbose = FALSE)
epi <- RunPCA(epi, npcs = 30, verbose = FALSE)
epi <- RunUMAP(epi, reduction = "pca", dims = 1:30)
epi <- RunTSNE(epi, reduction = "pca", dims = 1:30)
epi <- FindNeighbors(epi, dims = 1:30)
epi <- FindClusters(epi, resolution = 0.8)
DotPlot(epi, features = epi.markers, group.by = "integrated_snn_res.0.8", assay= "RNA") + RotatedAxis()
at1 <- c(0,6,8,9,10)
at2 <- c(1,2,3,5,7,11,13,14)
epi$cell.type <- ifelse(epi$integrated_snn_res.0.8 %in% at1, "AT1",
                        ifelse(epi$integrated_snn_res.0.8 %in% at2, "AT2", "Airway epithelial"))
DimPlot(epi, reduction = "umap", group.by = "cell.type")
subsce$cell.type[subsce$cell.type == "epi"] <- epi$cell.type


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
immu <- subset(subsce,cell.type == "imm")
immu <- ScaleData(immu, verbose = FALSE)
immu <- RunPCA(immu, npcs = 30, verbose = FALSE)
immu <- RunUMAP(immu, reduction = "pca", dims = 1:30)
immu <- RunTSNE(immu, reduction = "pca", dims = 1:30)
immu <- FindNeighbors(immu, dims = 1:30)
immu <- FindClusters(immu, resolution = 0.5)
DotPlot(immu, features = immu.markers, group.by = "integrated_snn_res.0.5",assay = "RNA") + RotatedAxis()
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(1,16,20)]<- "B cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(11,13,14)]<- "Plasma cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(0,2,3,7,10,15)] <- "T cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% 8] <- "NK cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(4,12,18)] <- "Neutrophil"
immu$cell.type[immu$integrated_snn_res.0.5 %in% 9] <- "Mast cell"
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(5,6)] <- "Macrophage"
immu$cell.type[immu$integrated_snn_res.0.5 %in% c(17,19)] <- "Dendritic cell"
DimPlot(immu, reduction = "tsne", group.by = "cell.type")
subsce$cell.type[subsce$cell.type == "imm"] <- immu$cell.type
sce$cell.type[sce$cell.type == "NULL"]<-subsce$cell.type

sp <- as.character(unique(sce$orig.ident))
for(i in sp){
  sample <- subset(sce, orig.ident == i)
  saveRDS(sample,paste(i,".rds",sep = ""))
}

allsample <- list.files(pattern =  ".rds$")[-25:-26]
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


newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "1"],split="_"),function(x){x[1]}),"_p1t1",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "2"],split="_"),function(x){x[1]}),"_p1t2",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "3"],split="_"),function(x){x[1]}),"_p1t3",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "4"],split="_"),function(x){x[1]}),"_p1t4",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "5"],split="_"),function(x){x[1]}),"_p2t1",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "6"],split="_"),function(x){x[1]}),"_p2t2",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "7"],split="_"),function(x){x[1]}),"_p3t1",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "8"],split="_"),function(x){x[1]}),"_p3t2",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "9"],split="_"),function(x){x[1]}),"_p4t1",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "10"],split="_"),function(x){x[1]}),"_p4t2",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "11"],split="_"),function(x){x[1]}),"_p4t3",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "12"],split="_"),function(x){x[1]}),"_p5t1",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "13"],split="_"),function(x){x[1]}),"_p5t2",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "14"],split="_"),function(x){x[1]}),"_p5t3",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "15"],split="_"),function(x){x[1]}),"_p6t1",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "16"],split="_"),function(x){x[1]}),"_p6t2",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "17"],split="_"),function(x){x[1]}),"_p7t1",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "18"],split="_"),function(x){x[1]}),"_p7t2",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "19"],split="_"),function(x){x[1]}),"_p8t1",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "20"],split="_"),function(x){x[1]}),"_p8t2",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "21"],split="_"),function(x){x[1]}),"_nsc035lt",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "22"],split="_"),function(x){x[1]}),"_nsc036lt",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "23"],split="_"),function(x){x[1]}),"_nsc037lt",sep=""))
newcol <- c(newcol,paste(lapply(strsplit(colnames(sce)[n == "24"],split="_"),function(x){x[1]}),"_nsc040lt",sep=""))
