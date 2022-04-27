run_NATMI_output = function(name){
  dir <- paste("natmi/",name,"/out/Network_exp_0_spe_0_det_0.2_top_0_signal_lrc2p_weight_mean/Edges.csv", sep = "")
  posi <- read.csv(dir)
  posi <- posi[,c(1,4,2,3)]
  colnames(posi) <- c("source","target","ligand","receptor")
  posi
}