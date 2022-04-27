CellPhoneDB_process = function(sample){
# save data
dir <- paste("cpdb/",substitute(sample),sep = "")
dir.create(dir)
writeMM(sample@assays$RNA@counts, file = paste(dir,"/matrix.mtx",sep=""))
write(x = rownames(sample@assays$RNA@counts), file = paste(dir,"/features.tsv",sep=""))
write(x = colnames(sample@assays$RNA@counts), file = paste(dir,"/barcodes.tsv",sep=""))

# save meta
sample$cell_type <- Idents(sample)
sample@meta.data$Cell = rownames(sample@meta.data)
df <- sample@meta.data[, c('Cell', 'cell_type')]
write.table(df, file = paste(dir,"/meta.tsv",sep=""), sep = '\t', quote = F, row.names = F)
}