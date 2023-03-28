sce <- CreateSeuratObject(counts = GSE135893_ILD_annotated_fullsize@assays$RNA@counts)
sce$orig.ident <- GSE135893_ILD_annotated_fullsize$orig.ident
sce$Diagnosis <- GSE135893_ILD_annotated_fullsize$Diagnosis
sce$Sample_Name <- GSE135893_ILD_annotated_fullsize$Sample_Name
sce$Status <- GSE135893_ILD_annotated_fullsize$Status
subsce <- subset(sce, Diagnosis == "Control" | Diagnosis == "IPF")
saveRDS(subsce, "sce.rds")
sce <- readRDS("sce.rds")
sce <- NormalizeData(sce, verbose = FALSE)
sce <- FindVariableFeatures(sce, selection.method = "vst")
sce <- ScaleData(sce, verbose = FALSE)
sce <- RunPCA(sce, npcs = 30, verbose = FALSE)
ElbowPlot(sce, ndims = 30)
sce<- RunUMAP(sce, reduction = "pca", dims = 1:30)
sce <- RunTSNE(sce, reduction = "pca", dims = 1:30)
sce <- FindNeighbors(sce, dims = 1:30)
sce <- FindClusters(sce, resolution = 0.5)
DimPlot(sce, reduction = "tsne", group.by = "seurat_clusters")

inn.markers <- c("EPCAM", # epithelium
                 "PECAM1", # Endothelial
                 "COL1A1","PDGFRA", # Fibroblast
                 "CNN1","ACTA2","TAGLN","RGS5","DES", # Smooth Muscle
                 "PTPRC","CD3G","CD3E","CD79A" # immune
                 )
DotPlot(sce, features = inn.markers, group.by = "RNA_snn_res.0.5") + RotatedAxis()
imm <- c("1","2","5","6","8","9","10","13","18","21","22","23","24")
epi <- c("0","3","4","7","15","17","25","26","28")
endo <- c("12","14","16","19")
fibro <- c("11","27")
sce$cell.type <- ifelse(sce$RNA_snn_res.0.5 %in% imm, "immune",
                    ifelse(sce$RNA_snn_res.0.5 %in% epi, "epi", 
                           ifelse(sce$RNA_snn_res.0.5 %in% endo, "Endothelial",
                                  ifelse(sce$RNA_snn_res.0.5 %in% fibro, "Fibroblast","Smooth muscle"))))
DimPlot(sce, reduction = "umap", group.by = "cell.type")



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
at1 <- c("3","5","9")
at2 <- c("1","2","4","6","7","8","10")
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
immu$cell.type[immu$RNA_snn_res.0.5 %in% 18] <- "B cell"
immu$cell.type[immu$RNA_snn_res.0.5 %in% c(20,23)] <- "Plasma cell"
immu$cell.type[immu$RNA_snn_res.0.5 %in% c(5,8,16,21)] <- "T cell"
immu$cell.type[immu$RNA_snn_res.0.5 %in% 11] <- "NK cell"
immu$cell.type[immu$RNA_snn_res.0.5 %in% c(10,13)] <- "Neutrophil"
immu$cell.type[immu$RNA_snn_res.0.5 %in% 19] <- "Mast cell"
immu$cell.type[immu$RNA_snn_res.0.5 %in% c(1,2,4,9,14,15)] <- "Macrophage"
immu$cell.type[immu$RNA_snn_res.0.5 %in% c(3,22)] <- "Dendritic cell"
immu$cell.type[immu$RNA_snn_res.0.5 %in% c(0,6,7,12,17)] <- "Monocyte"
DimPlot(immu, reduction = "umap", group.by = "cell.type")
sce$cell.type[sce$cell.type == "immune"] <- immu$cell.type

ctlsce <- subset(sce, Diagnosis == "Control")
ipfsce <- subset(sce, Diagnosis == "IPF")

ctlsp <- unique(ctlsce$Sample_Name)
for(i in ctlsp){
  sp <- subset(ctlsce, Sample_Name == i)
  saveRDS(sp,paste("healthy/",i,".rds",sep = ""))
}
ipfsp <- unique(ipfsce$Sample_Name)
for(i in ipfsp){
  sp <- subset(ipfsce, Sample_Name == i)
  saveRDS(sp,paste("ipf/",i,".rds",sep = ""))
}
