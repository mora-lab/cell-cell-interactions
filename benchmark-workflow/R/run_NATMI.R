run_NATMI_output = function(name,inputdata){
  dir <- paste("Python/natmi/",name,"/out/Network_exp_0_spe_0_det_0.2_top_0_signal_lrc2p_weight_mean/Edges.csv", sep = "")
  time1 <- read.table(paste("Python/natmi/",name,"/out/time1.txt", sep = ""), row.names = 1)
  time2 <- read.table(paste("Python/natmi/",name,"/out/time2.txt", sep = ""), row.names = 1)
  speed <- data.frame(cells=ncol(inputdata), time=time1["real",1]+time2["real",1], row.names = NULL)
  result <- suppressWarnings(try(read.csv(dir), silent = TRUE))
  if('try-error' %in% class(result)){
    posi <- data.frame(source=NA, target=NA, ligand=NA, receptor=NA)
    list(lrpairs=posi, speed=speed, pairs=0)
  }else{posi <- read.csv(dir)
    posi <- posi[,c(1,4,2,3)]
    colnames(posi) <- c("source","target","ligand","receptor")
    list(lrpairs=posi, speed=speed, pairs=length(rownames(posi)))
  }
}
