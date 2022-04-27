NATMI_process = function(sample){
  # save data
  dir <- paste("natmi/",substitute(sample),sep = "")
  dir.create(dir)
  write.csv(as.matrix(GetAssayData(object = sample, assay = "RNA", slot = "data")), paste(dir,"/em.csv",sep=""), row.names = T)
  # save meta
  write.csv(Idents(sample), paste(dir,"/ann.csv",sep=""), row.names = T)
}
