run_scMLnet_SSP = function(posi,sample,gold){
database <- read.table("https://raw.githubusercontent.com/mora-lab/cell-cell-interactions/main/benchmark-workflow/R/scMLnet_database/LigRec.txt", header = TRUE, sep = "\t")
database <- database[database$Ligand %in% gold$ligand,]
database <- database[database$Receptor %in% gold$receptor,]
# create database with all cell pairs
types <- as.character(unique(Idents(sample)))
singlelist <- data.frame(source=NA, target=NA)
full.database <- data.frame()

for(i in 1:length(types)){
  cellA <- types[i]
  for(n in 1:length(types)){
    cellB <- types[n]
    singlelist[c(1:nrow(database)),1] <- cellA
    singlelist[c(1:nrow(database)),2] <- cellB
    singlelist <- cbind(singlelist, database[,1])
    full.database <- rbind(full.database, singlelist)
    singlelist <- data.frame(source=NA, target=NA)
  }
}

colnames(full.database)[3] <- "ligrec"

# create all false results list
posi[,4] <- paste(posi[,1],posi[,2],posi[,3],sep = "-")
full.database[,4] <- paste(full.database[,1],full.database[,2],full.database[,3],sep = "-")
falselist <- full.database[!full.database$V4 %in% posi$V4,]

# create gold standard list
gold[,5] <- paste(gold[,1],gold[,2],gold[,3],sep = "-")
gold[,6] <- paste(gold[,5],gold[,4],sep = "_")

# calculate TP,FP,TN,FN
TP <- nrow(posi[posi$V4 %in% gold$V6,])
FP <- nrow(posi) - TP
TN <- nrow(falselist[!falselist$V4 %in% gold$V6,])
FN <- nrow(falselist) - TN

# calculate accuracy, sensitivity and specificity
acc <- TP/(TP+FP)
sst <- TP/(TP+FN)
spc <- TN/(FP+TN)
scMLnetSSP <- list(Accuracy=acc, Sensitivity=sst, Specificity=spc)
scMLnetSSP
}