run_CellChat = function(inputdata){
# preparation of input data
data <- inputdata@assays$RNA@data
meta <- data.frame(inputdata@active.ident)
f = file()
sink(file = f, type = c("output","message"))
s <- proc.time()
CellChat <- createCellChat(object = data, meta = meta, group.by = colnames(meta))
# choose database
CellChatDB <- CellChatDB.human
CellChat@DB <- CellChatDB

# processing data
CellChat <- subsetData(CellChat)
CellChat <- identifyOverExpressedGenes(CellChat)
CellChat <- identifyOverExpressedInteractions(CellChat)

# calculate the network
CellChat <- computeCommunProb(CellChat)
result <- try(subsetCommunication(CellChat), silent = TRUE)
e <- proc.time()
speed <- data.frame(cells=ncol(inputdata), time=(e-s)[3], row.names = NULL)
sink()
close(f)
# convert into unified format
if('try-error' %in% class(result)){
  full.cclist <- data.frame(source=NA, target=NA, ligand=NA, receptor=NA)
  list(lrpairs=full.cclist, pairs=0, speed=speed)
}else{
 cclist <- result[,c(1,2,3,4,5,6,7)]
 cclist[,7] <- as.character(cclist[,7])
 full.cclist <- data.frame()
 for (i in 1:nrow(cclist)) {
  if("_" %in% strsplit(cclist[i,4],split = "")[[1]] || 
     "complex" %in% strsplit(cclist[i,4],split = " ")[[1]] ||
     "receptor" %in% strsplit(cclist[i,4],split = " ")[[1]]){
    lr <- strsplit(cclist[i,7],split = "_")[[1]]
    for (j in 2:length(lr)) {
      newlr <- cclist[i,]
      newlr[1,4] <- lr[j]
      full.cclist <- rbind(full.cclist, newlr)
    }
  }
  else{ full.cclist <- rbind(full.cclist, cclist[i,])
  }
}
full.cclist <- unique(full.cclist[,c(1:6)])
full.cclist[,1] <- as.character(full.cclist[,1])
full.cclist[,2] <- as.character(full.cclist[,2])
list(lrpairs=full.cclist, pairs=nrow(full.cclist), speed=speed)
}
}
