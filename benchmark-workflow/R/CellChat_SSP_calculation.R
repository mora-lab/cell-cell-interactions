# create database without complex pairs
# CellChatDB <- CellChatDB.human
# dblist <- CellChatDB$interaction[,c(3,4)]
# full.dblist <- data.frame()
# for (i in 1:nrow(dblist)) {
#     if("_" %in% strsplit(dblist[i,2],split = "")[[1]]){
#        lr <- strsplit(row.names(dblist[i,]),split = "_")[[1]]
#        for (j in 2:length(lr)) {
#           newlr <- dblist[i,]
#           newlr[1,2] <- lr[j]
#           full.dblist <- rbind(full.dblist, newlr)
#        }
#    }
#     else{ full.dblist <- rbind(full.dblist, dblist[i,])
#    }
# }
# modify rest complex pairs manually

run_CellChat_SSP = function(full.cclist, gold){
full.cclist <- full.cclist[,1:4]
full.dblist <- read.csv("git@github.com:mora-lab/cell-cell-interactions/blob/main/benchmark-workflow/R/CellChatDatabase.csv", header = TRUE, row.names = 1)
# create database with all cell pairs
rownames(full.dblist) <- NULL
types <- as.character(unique(Idents(sample)))
singlelist <- data.frame(source=NA, target=NA)
full.database <- data.frame()

for(i in 1:length(types)){
cellA <- types[i]
 for(n in 1:length(types)){
   cellB <- types[n]
   singlelist[c(1:nrow(full.dblist)),1] <- cellA
   singlelist[c(1:nrow(full.dblist)),2] <- cellB
   singlelist <- cbind(singlelist, full.dblist)
   full.database <- rbind(full.database, singlelist)
   singlelist <- data.frame(source=NA, target=NA)
 }
}

# create all false results list
full.cclist[,5] <- paste(full.cclist[,1],full.cclist[,2],full.cclist[,3],full.cclist[,4],sep = "-")
full.database[,5] <- paste(full.database[,1],full.database[,2],full.database[,3],full.database[,4],sep = "-")
falselist <- full.database[!full.database$V5 %in% full.cclist$V5,]

# create gold standard list
gold[,5] <- paste(gold[,1],gold[,2],gold[,3],gold[,4],sep = "-")

# calculate TP,FP,TN,FN
TP <- nrow(full.cclist[full.cclist$V5 %in% gold$V5,])
FP <- nrow(full.cclist) - TP
TN <- nrow(falselist[!falselist$V5 %in% gold$V5,])
FN <- nrow(falselist) - TN

# calculate accuracy, sensitivity and specificity
acc <- TP/(TP+FP)
sst <- TP/(TP+FN)
spc <- TN/(FP+TN)
CellChatSSP <- list(Accuracy=acc, Sensitivity=sst, Specificity=spc)
CellChatSSP
}


