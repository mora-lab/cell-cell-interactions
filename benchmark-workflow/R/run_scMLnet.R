run_scMLnet = function(inputdata){
mtx <- inputdata@assays$RNA@counts
inputdata$Cluster <- as.character(Idents(inputdata))
inputdata@meta.data$Barcode = rownames(inputdata@meta.data)
ann <- inputdata@meta.data[, c('Barcode', 'Cluster')]
write.table(ann, file ='meta.tsv', sep = '\t', quote = F, row.names = F)
annfile <- "meta.tsv"
types <- unique(ann$Cluster)
maxtype <- length(types)
pval <- 0.05
logfc <- 0.15
LigRecLib <- "git@github.com:mora-lab/cell-cell-interactions/blob/main/benchmark-workflow/R/scMLnet_database/LigRec.txt"
TFTarLib <- "git@github.com:mora-lab/cell-cell-interactions/blob/main/benchmark-workflow/R/scMLnet_database/TFTargetGene.txt"
RecTFLib <- "git@github.com:mora-lab/cell-cell-interactions/blob/main/benchmark-workflow/R/scMLnet_database/RecTF.txt"
netList <- NULL
for (i in 1:maxtype) {
LigClu <- types[i]
    for(l in 1:maxtype) {
      if (l == i) {}
      else{
      RecClu <- types[l]
      net <- try(RunMLnet(mtx, annfile, RecClu, LigClu, pval, logfc, LigRecLib, TFTarLib, RecTFLib), silent=FALSE)
      if('try-error' %in% class(net))           
      {net <- NA}
      net <- list(net)
      names(net) <- paste(LigClu,"-",RecClu)
      netList <- c(netList,net)
      }
   }
}
posi <- data.frame()
singlelist <- data.frame(source=NA, target=NA, ligrec=NA)
for (i in 1:length(netList)) {
   if (!netList[[i]] %in% NA) {
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
posi
}