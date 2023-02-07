run_SSP = function(posi,allparis,gold){
  posi <- posi[,1:4]
  posi[,5] <- paste(posi[,1],posi[,2],posi[,3],posi[,4],sep = "-")
  gold[,5] <- paste(gold[,1],gold[,2],gold[,3],gold[,4],sep = "-")
  
  # calculate TP,FP,TN,FN
  TP <- as.numeric(nrow(posi[posi$V5 %in% gold$V5,]))
  FP <- as.numeric(nrow(posi) - TP)
  FN <- as.numeric(nrow(gold) - TP)
  TN <- as.numeric(nrow(allpairs) - (TP + FP + FN))
  
  # calculate accuracy, sensitivity and specificity
  acc <- (TP+TN)/nrow(allpairs)
  pre <- TP/(TP+FP)
  sst <- TP/(TP+FN)
  spc <- TN/(FP+TN)
  f1 <- 2*TP/(2*TP+FP+FN)
  mcc <- (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  npv <- TN/(TN+FN)
  SSP <- list(Precision=pre, Sensitivity=sst, Specificity=spc, F1score=f1, MCC=mcc, NPV=npv)
  SSP
}

