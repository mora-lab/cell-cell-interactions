library(Seurat)
library(stringr)
library(dplyr)

sce <- readRDS("GSE135893_ILD_annotated_fullsize.rds")
levels(Idents(sce))
new.cluster.ids <- c("AT1","AT2","B Cells","Basal","Dendritic cell","Ciliated","Differentiating Ciliated",
                     "Endothelial","Fibroblast","HAS1 High Fibroblasts","KRT5-/KRT17+",
                     "Lymphatic Endothelial Cells","Macrophage","Mast","Mesothelial Cells",             
                     "Monocyte","MUC5AC+ High","MUC5B+","Myofibroblasts","NK Cells","Dendritic cell",                          
                     "Plasma Cells","PLIN2+ Fibroblasts","Proliferating Epithelial Cells",
                     "Proliferating Macrophages","Proliferating T Cells","SCGB3A2+",                      
                     "SCGB3A2+ SCGB1A1+","Smooth Muscle Cells","Tcell","Transitional AT2")
names(new.cluster.ids) <- levels(sce)
sce <- RenameIdents(sce, new.cluster.ids)
table(Idents(sce))
sce$cell.type <- as.character(Idents(sce))
sce$Class <- sce$Diagnosis

# Simplify data
IPFsample <- subset(sce, Class == "IPF")
cell <- colnames(IPFsample[,Idents(IPFsample) %in% c("Monocyte","Macrophage","Fibroblast","Endothelial","Tcell","AT2","AT1","Mast")])
IPFsample <- subset(IPFsample, cells = cell) 
IPFsample2 <- CreateSeuratObject(counts = IPFsample@assays$RNA, meta.data = IPFsample@meta.data)
IPFsample2 <- NormalizeData(IPFsample2, verbose = FALSE)
Idents(IPFsample2) <- "cell.type"
sampledata <- unique(IPFsample2$Sample_Name)
for(x in sampledata){
  sub <- subset(IPFsample2, Sample_Name == x)
  saveRDS(sub,paste("GSE135893/", x, ".rds", sep = ""))
}
