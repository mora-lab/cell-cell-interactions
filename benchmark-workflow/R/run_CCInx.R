run_CCInx = function(inputdata,top){
gsl <- BuildGeneStatList(inD = inputdata, cl = inputdata@active.ident, assayType = "RNA")
inx <- BuildCCInx(GeneStatList=gsl)
posi <- inx$edges
print(paste("Total predicted L-R pairs:",nrow(posi)))
posi <- posi[order(posi$edgeWeight, decreasing=T),][1:top,]
for (i in 1:nrow(posi)) {
  lig <- strsplit(posi[i,1],split = "_")[[1]]
  rec <- strsplit(posi[i,2],split = "_")[[1]]
  posi[i,1:5] <- c(lig[2], rec[2], lig[1], rec[1], posi[i,3])
}
colnames(posi) <- c("source","target","ligand","receptor", "edgeWeight")
rownames(posi) <- NULL

print(paste("Choose top L-R pairs:",top))
posi
}