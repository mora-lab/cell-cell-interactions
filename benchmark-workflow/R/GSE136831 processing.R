library(Seurat)
library(stringr)
library(dplyr)

matrix <- ReadMtx("matrix.mtx.gz","barcodes.txt.gz","features.txt.gz",skip.feature = 1)
meta <- read.table("GSE136831_AllCells.Samples.CellType.MetadataTable.txt", header = T, sep = "\t")
sce <- CreateSeuratObject(counts = matrix)
sce$Class <- meta$Disease_Identity
sce$orig.ident <- meta$Subject_Identity
sce$cell.type <- meta$Manuscript_Identity
sce <- NormalizeData(sce, verbose = FALSE)
Idents(sce) <- "cell.type"
levels(Idents(sce))
new.cluster.ids <- c("Monocyte","Alveolar macrophage","NK cell","Monocyte","Lymphatic","Macrophage",
                     "Fibroblast","Basal cell","Myofibroblast","Multiplet","Dendritic cell","B cell",
                     "Endothelial","Tcell","AT2","AT1","Endothelial","Endothelial",
                     "Ciliated cell","Club cell","Endothelial","Plasma","ILC","T cytotoxic",
                     "Goblet cell","T regulatory","Mesothelial","ILC","SMC","Dendritic cell",
                     "Dendritic cell","Mast","Dendritic cell","Dendritic cell","Endothelial","Pericyte",
                     "Aberrant basaloid","Ionocyte","PNEC")
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)
table(Idents(sce))
sce$cell.type <- as.character(Idents(sce))

IPFsample <- subset(sce, Class == "IPF")
cell <- colnames(IPFsample[,Idents(IPFsample) %in% c("Monocyte","Macrophage","Fibroblast","Endothelial","Tcell","AT2","AT1","Mast")])
IPFsample <- subset(IPFsample, cells = cell) 
sampledata <- unique(IPFsample$orig.ident)
for(x in sampledata){
  sub <- subset(IPFsample, ori.ident == x)
  saveRDS(sub,paste("GSE136831/", x, ".rds", sep = ""))
}