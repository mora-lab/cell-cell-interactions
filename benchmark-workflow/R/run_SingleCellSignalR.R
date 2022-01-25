run_SingleCellSignalR = function(inputdata){
data <- data.frame(inputdata@assays$RNA@data, stringsAsFactors = F)
cluster <- as.integer(Idents(inputdata))
c.names <- levels(Idents(inputdata))
signal <- cell_signaling(data = data, genes = rownames(data), cluster = cluster, c.names = c.names, write = FALSE)
posi <- data.frame()
singlelist <- data.frame(source=NA, target=NA, ligand=NA, receptor=NA, LRscore=NA)
for (i in 1:length(signal)) {
  pair <- data.frame(signal[i])[,c(1,2,4)]
  CellA <- strsplit(names(signal)[i],split = "-")[[1]][1]
  CellB <- strsplit(names(signal)[i],split = "-")[[1]][2]
  for (j in 1:nrow(pair)) {
    singlelist[j,] <- c(CellA, CellB, pair[j,1],pair[j,2], pair[j,3])
  }
  posi <- rbind(posi, singlelist)
  singlelist <- data.frame(source=NA, target=NA, ligand=NA, receptor=NA, LRscore=NA)
}
posi
}