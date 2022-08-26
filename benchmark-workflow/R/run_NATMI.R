run_NATMI_output = function(name){
  dir <- paste("Python/natmi/",name,"/out/Network_exp_0_spe_0_det_0.2_top_0_signal_lrc2p_weight_mean/Edges.csv", sep = "")
  result <- suppressWarnings(try(read.csv(dir), silent = TRUE))
  if('try-error' %in% class(result)){
    print("Total predicted L-R pairs: 0")
    posi <- data.frame(source=NA, target=NA, ligand=NA, receptor=NA)
    posi
  }else{posi <- read.csv(dir)
    posi <- posi[,c(1,4,2,3)]
    colnames(posi) <- c("source","target","ligand","receptor")
    print(paste("Total predicted L-R pairs:",nrow(posi)))
    posi
  }
}