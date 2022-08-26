CellPhoneDB_process = function(seu,sample){
# save data
dir <- paste("Python/cpdb/", sample, sep = "")
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