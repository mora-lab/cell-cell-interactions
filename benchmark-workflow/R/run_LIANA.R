run_Liana = function(inputdata,top){
  options(warn = -1)
  f = file()
  sink(file = f, type = c("output","message"))
  s <- proc.time()
  liana <- liana_wrap(inputdata,assay = "RNA")
  liana <- liana_aggregate(liana)
  e <- proc.time()
  speed <- data.frame(cells=ncol(inputdata), time=(e-s)[3], row.names = NULL)
  sink()
  close(f)
  liana <- as.data.frame(liana[order(liana$aggregate_rank,decreasing=FALSE),1:5][1:top,])
  posi <- data.frame()
  for(i in 1:top){
    l <- strsplit(liana[i,3],split="_")[[1]]
    r <- strsplit(liana[i,4],split="_")[[1]]
    lr <- expand.grid(l,r)
    st <- data.frame(source=rep(liana[i,1],time=nrow(lr)),target=rep(liana[i,2],time=nrow(lr)))
    st <- cbind(st,lr)
    posi <- rbind(posi,st)
  }
  colnames(posi)[3:4] <- c("ligand","receptor")
  posi[,3] <- as.character(posi[,3])
  posi[,4] <- as.character(posi[,4])
  list(lrpairs=posi, speed=speed, pairs=nrow(posi))
}