NATMI_process = function(seu, sample){
  # save data
  dir <- paste("Python/natmi/", sample,sep = "")
  dir.create(dir)
  write.csv(as.matrix(GetAssayData(object = seu, assay = "RNA", slot = "data")), paste(dir,"/em.csv",sep=""), row.names = T)
  # save meta
  write.csv(Idents(seu), paste(dir,"/ann.csv",sep=""), row.names = T)
}
