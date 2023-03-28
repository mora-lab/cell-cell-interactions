suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Matrix))
allsample <- c(paste("healthy/",list.files("healthy/"),sep=""), 
               paste("ipf/",list.files("ipf/"),sep=""))
dir.create("healthy/cpdb")
for (sample in allsample){
  seu <- readRDS(sample)
  Idents(seu) <- "cell.type"
  CellPhoneDB_process(seu,"healthy/cpdb/")
}

CellPhoneDB_process = function(seu,dir){
  # save data
  sample <- unique(seu$Sample_Name)
  dir <- paste(dir, sample, sep = "")
  dir.create(dir)
  writeMM(seu@assays$RNA@counts, file = paste(dir,"/matrix.mtx",sep=""))
  write(x = rownames(seu@assays$RNA@counts), file = paste(dir,"/features.tsv",sep=""))
  write(x = colnames(seu@assays$RNA@counts), file = paste(dir,"/barcodes.tsv",sep=""))
  
  # save meta
  seu$cell_type <- Idents(seu)
  seu@meta.data$Cell = rownames(seu@meta.data)
  df <- seu@meta.data[, c('Cell', 'cell_type')]
  write.table(df, file = paste(dir,"/meta.tsv",sep=""), sep = '\t', quote = F, row.names = F)
}
