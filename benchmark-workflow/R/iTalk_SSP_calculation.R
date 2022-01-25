run_iTalk_SSP = function(posi,gold){
load("git@github.com:mora-lab/cell-cell-interactions/blob/main/benchmark-workflow/R/iTalk_LR_database.rda")
database <- database[,c(2,4)]
colnames(database) <- c("ligand","receptor")
posi <- posi[,1:4]

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
    singlelist <- cbind(singlelist, database)
    full.database <- rbind(full.database, singlelist)
    singlelist <- data.frame(source=NA, target=NA)
  }
}

# create all false results list
posi[,5] <- paste(posi[,1],posi[,2],posi[,3],posi[,4],sep = "-")
full.database[,5] <- paste(full.database[,1],full.database[,2],full.database[,3],full.database[,4],sep = "-")
falselist <- full.database[!full.database$V5 %in% posi$V5,]

# create gold standard list
gold[,5] <- paste(gold[,1],gold[,2],gold[,3],gold[,4],sep = "-")

# calculate TP,FP,TN,FN
TP <- nrow(posi[posi$V5 %in% gold$V5,])
FP <- nrow(posi) - TP
TN <- nrow(falselist[!falselist$V5 %in% gold$V5,])
FN <- nrow(falselist) - TN

# calculate accuracy, sensitivity and specificity
acc <- TP/(TP+FP)
sst <- TP/(TP+FN)
spc <- TN/(FP+TN)
iTalkSSP <- list(Accuracy=acc, Sensitivity=sst, Specificity=spc)
iTalkSSP
}