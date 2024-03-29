run_iTalk = function(inputdata,top){
# data preparation
data <- data.frame(inputdata@assays$RNA@data, stringsAsFactors = F)
cell_type <- as.character(inputdata@active.ident)
data <- t(data)
data <- as.data.frame(data)
data <- cbind(data, cell_type)
s <- proc.time()
# find top 50 percent highly expressed genes
highly_exprs_genes <-rawParse(data, top_genes=50, stats='mean')
# find the ligand-receptor pairs from highly expressed genes
comm_list <-c('growth factor','other','cytokine','checkpoint')
res <- NULL

for(comm_type in comm_list){
  res_cat <- try(FindLR( highly_exprs_genes, datatype='mean count', comm_type=comm_type), silent=TRUE)
  res_cat <- res_cat[ order( res_cat$cell_from_mean_exprs*res_cat$cell_to_mean_exprs, decreasing=T),]
  res <-rbind(res,res_cat)
}
e <- proc.time()
rm(data)
gc(verbose=FALSE)
speed <- data.frame(cells=ncol(inputdata), time=(e-s)[3], row.names = NULL)
if(nrow(res)>0){
  posi <- res[,c(4,6,1,2,3,5)]
  posi <- posi[order(posi$cell_from_mean_exprs*posi$cell_to_mean_exprs,decreasing=T),][1:top,]
  colnames(posi) <- c("source","target","ligand","receptor","source_mean_exprs","target_mean_exprs")
  list(lrpairs=posi, speed=speed, pairs=nrow(posi))
}else{
  posi <- data.frame(source=NA, target=NA, ligand=NA, receptor=NA,source_mean_exprs=NA,target_mean_exprs=NA)
  list(lrpairs=posi, speed=speed, pairs=0)
}
}