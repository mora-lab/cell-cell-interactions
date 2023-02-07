run_CCInx = function(inputdata,top){
s <- proc.time()
f = file()
sink(file = f, type = c("output","message"))
gsl <- try(BuildGeneStatList(inD = inputdata, cl = inputdata@active.ident, assayType = "RNA"), silent = TRUE)
if('try-error' %in% class(gsl)){
sink()
close(f)
posi <- data.frame(source=NA, target=NA, ligand=NA, receptor=NA, edgeWeight=NA)
speed <- data.frame(cells=ncol(inputdata), time=NA)
list(lrpairs=posi, speed=speed, pairs=0)
}else{
inx <- BuildCCInx(GeneStatList=gsl)
e <- proc.time()
speed <- data.frame(cells=ncol(inputdata), time=(e-s)[3], row.names = NULL)
sink()
close(f)
ligand <- inx$nodes[inx$nodes$proteinType %in% c("Ligand","ECM/Ligand","ECM/Receptor/Ligand","Receptor/Ligand"),"node"]
receptor <- inx$nodes[inx$nodes$proteinType %in% c("Receptor","ECM/Receptor","ECM/Receptor/Ligand","Receptor/Ligand"),"node"]
lr <- inx$edges[(inx$edges$nodeA %in% ligand) & (inx$edges$nodeB %in% receptor),]
rl <- inx$edges[(inx$edges$nodeB %in% ligand) & (inx$edges$nodeA %in% receptor),c("nodeB","nodeA","edgeWeight")]
colnames(rl) <- c("nodeA","nodeB","edgeWeight")
posi <- rbind(lr,rl)
posi <- posi[order(posi$edgeWeight, decreasing=T),][1:top,]
for (i in 1:nrow(posi)) {
  lig <- strsplit(posi[i,1],split = "_")[[1]]
  rec <- strsplit(posi[i,2],split = "_")[[1]]
  posi[i,1:5] <- c(lig[2], rec[2], lig[1], rec[1], posi[i,3])
}
colnames(posi) <- c("source","target","ligand","receptor", "edgeWeight")
rownames(posi) <- NULL
list(lrpairs=posi, speed=speed, pairs=nrow(posi))
}
}