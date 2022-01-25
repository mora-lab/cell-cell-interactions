run_iTalk = function(inputdata,threshold){
# data preparation
data <- data.frame(inputdata@assays$RNA@data, stringsAsFactors = F)
cell_type <- as.character(inputdata@active.ident)
data <- t(data)
data <- as.data.frame(data)
data <- cbind(data, cell_type)

# find top 50 percent highly expressed genes
highly_exprs_genes <-rawParse(data, top_genes=50, stats='mean')

# find the ligand-receptor pairs from highly expressed genes
comm_list <-c('growth factor','other','cytokine','checkpoint')
res <- NULL

for(comm_type in comm_list){
  res_cat <- FindLR( highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  res_cat <- res_cat[ order( res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs, decreasing=T),]
  res <-rbind(res,res_cat)
}
posi <- res[,c(4,6,1,2,3,5)]
posi <- posi[order(posi$cell_from_mean_exprs*posi$cell_to_mean_exprs,decreasing=T),][posi$cell_from_mean_exprs*posi$cell_to_mean_exprs > threshold,] # how to define the threshold? 
colnames(posi) <- c("source","target","ligand","receptor","source_mean_exprs","target_mean_exprs")
posi
}