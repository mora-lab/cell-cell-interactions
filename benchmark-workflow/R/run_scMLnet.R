run_scMLnet = function(inputdata){
options(warn = -1)
mtx <- inputdata@assays$RNA@counts
inputdata$Cluster <- as.character(Idents(inputdata))
inputdata@meta.data$Barcode = rownames(inputdata@meta.data)
ann <- inputdata@meta.data[, c('Barcode', 'Cluster')]
write.table(ann, file ='R/meta.tsv', sep = '\t', quote = F, row.names = F)
annfile <- "R/meta.tsv"
types <- unique(ann$Cluster)
maxtype <- length(types)
pval <- 0.05
logfc <- 0.15
LigRecLib <- "https://raw.githubusercontent.com/mora-lab/cell-cell-interactions/main/benchmark-workflow/R/scMLnet_database/LigRec.txt"
TFTarLib <- "https://raw.githubusercontent.com/mora-lab/cell-cell-interactions/main/benchmark-workflow/R/scMLnet_database/TFTargetGene.txt"
RecTFLib <- "https://raw.githubusercontent.com/mora-lab/cell-cell-interactions/main/benchmark-workflow/R/scMLnet_database/RecTF.txt"
netList <- NULL
f = file()
sink(file = f, type = c("output","message"))
s <- proc.time()
for (i in 1:maxtype) {
LigClu <- types[i]
    for(l in 1:maxtype) {
      if (l == i) {}
      else{
      RecClu <- types[l]
      net <- try(RunMLnet(mtx, annfile, RecClu, LigClu, pval, logfc, LigRecLib, TFTarLib, RecTFLib), silent = TRUE)
      if('try-error' %in% class(net))           
      {net <- NA}
      net <- list(net)
      names(net) <- paste(LigClu,"-",RecClu)
      netList <- c(netList,net)
      }
   }
}
e <- proc.time()
sink()
close(f)
speed <- data.frame(cells=ncol(inputdata), time=(e-s)[3], row.names = NULL)
posi <- data.frame()
singlelist <- data.frame(source=NA, target=NA, ligrec=NA)
for (i in 1:length(netList)) {
   if (!(netList[[i]][1] %in% NA | is.null(netList[[i]][[1]]))) {
      pair <- netList[[i]][[1]]
      CellA <- strsplit(names(netList)[i],split = " - ")[[1]][1]
      CellB <- strsplit(names(netList)[i],split = " - ")[[1]][2]
      for (j in 1:length(pair)) {
         singlelist[j,] <- c(CellA, CellB, pair[j])
      }
      posi <- rbind(posi, singlelist)
      singlelist <- data.frame(source=NA, target=NA, ligrec=NA)
   }
}

if(nrow(posi) == 0){
posi <- data.frame(source=NA, target=NA, ligand=NA, receptor= NA)
list(lrpairs=posi, pairs=0, speed=speed)
}else{
  lr <- strsplit(posi$ligrec, split = "_")
  colnames(posi) <- c("source", "target", "ligand")
  for(i in 1:nrow(posi)){
    posi[i,"ligand"] <- lr[[i]][1]
    posi[i,"receptor"] <- lr[[i]][2]
  }
  list(lrpairs=posi, pairs=nrow(posi), speed=speed)}
}