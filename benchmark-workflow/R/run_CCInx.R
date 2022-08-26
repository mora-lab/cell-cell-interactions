run_CCInx = function(inputdata,top){
predict <- try(BuildGeneStatList(inD = inputdata, cl = inputdata@active.ident, assayType = "RNA"), silent = TRUE)
if('try-error' %in% class(predict)){
print("Total predicted L-R pairs: 0")
posi <- data.frame(source=NA, target=NA, ligand=NA, receptor=NA, edgeWeight=NA)
posi
}else{s <- proc.time()
f = file()
sink(file = f, type = c("output","message"))
gsl <- BuildGeneStatList(inD = inputdata, cl = inputdata@active.ident, assayType = "RNA")
inx <- BuildCCInx(GeneStatList=gsl)
e <- proc.time()
sink()
close(f)
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
print("Time comsuming:")
print(e-s)
posi
}
}