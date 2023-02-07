run_CellPhoneDB_output = function(name,inputdata){
dir <- paste("Python/cpdb/",name,"/out/significant_means.txt", sep = "")
posi <- read.table(dir, header = T, sep = "\t")
time <- read.table(paste("Python/cpdb/",name,"/out/time.txt", sep = ""), row.names = 1)
speed <- data.frame(cells=ncol(inputdata), time=time["real",1], row.names = NULL)
list <- data.frame()

for (i in 13:ncol(posi)) {
 pval <- posi[!posi[,i] %in% NA,c(3:6,i)]
 if(nrow(pval)>0){
   cells <- strsplit(colnames(pval)[5],split = "\\.")[[1]]
   for (n in 1:nrow(pval)) {
   new <- as.character(c(pval[n,1:4], cells[1], cells[2]))
   list <- rbind(list, new)
   }
 }
}
colnames(list) <- c("partner_a","partner_b","gene_a","gene_b","cell_a","cell_b")

singleinter <- list[!(list$gene_a == ""|list$gene_b == ""),]
output <- NULL
for (i in 1:nrow(singleinter)) {
  output$cellA[i] <- singleinter[i,5]
  output$cellB[i] <- singleinter[i,6]
  output$geneA[i] <- singleinter[i,3]
  output$geneB[i] <- singleinter[i,4]
}
output <- as.data.frame(output)

complexinter <- list[list$gene_a == ""|list$gene_b == "",]
complex <- read.csv("https://raw.githubusercontent.com/mora-lab/cell-cell-interactions/main/benchmark-workflow/R/CellPhoneDB database/complex_input.csv")
gene <- read.csv("https://raw.githubusercontent.com/mora-lab/cell-cell-interactions/main/benchmark-workflow/R/CellPhoneDB database/gene_input.csv")
# to simplify codes, some changes were applied to CellPhoneDB's database and prediction results:
complex[17,1] <- "CD94%NKG2A"
complex[18,1] <- "CD94%NKG2C"
complex[19,1] <- "CD94%NKG2E"
try(complexinter[grep("^complex:CD94",complexinter$partner_a),1] <- paste("complex:CD94%", strsplit(complexinter[grep("^complex:CD94",complexinter$partner_a),1],split=":")[[1]][3],sep = ""),silent = TRUE)
try(complexinter[grep("^complex:CD94",complexinter$partner_b),2] <- paste("complex:CD94%", strsplit(complexinter[grep("^complex:CD94",complexinter$partner_b),2],split=":")[[1]][3],sep = ""),silent = TRUE)

for (i in 1:nrow(complexinter)) {
  if (complexinter[i,3] == "" & complexinter[i,4] == "") {
    proA <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],2],3][1]
    proB <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],3],3][1]
    proD <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],2],3][1]
    proE <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],3],3][1]
    output <- rbind(output,c(complexinter$cell_a[i],complexinter$cell_b[i],proA,proD),c(complexinter$cell_a[i],complexinter$cell_b[i],proA,proE),c(complexinter$cell_a[i],complexinter$cell_b[i],proB,proD),c(complexinter$cell_a[i],complexinter$cell_b[i],proB,proE))
    if (!complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],4] == ""){
      proC <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],4],3][1]
      output <- rbind(output,c(complexinter$cell_a[i],complexinter$cell_b[i],proC,proD),c(complexinter$cell_a[i],complexinter$cell_b[i],proC,proE))}
    if (!complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],4] == ""){
      proF <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],4],3][1]
      output <- rbind(output,c(complexinter$cell_a[i],complexinter$cell_b[i],proA,proF),c(complexinter$cell_a[i],complexinter$cell_b[i],proB,proF))}
    if (!complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],4] == "" & !complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],4] == ""){
      output <- rbind(output,c(complexinter$cell_a[i],complexinter$cell_b[i],proC,proF))}
  }else{
    if (complexinter[i,3] == "") {
      proA <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],2],3][1]
      proB <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],3],3][1]
      proD <- complexinter$gene_b[i]
      output <- rbind(output,c(complexinter$cell_a[i],complexinter$cell_b[i],proA,proD),c(complexinter$cell_a[i],complexinter$cell_b[i],proB,proD))
      if (!complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],4] == ""){
        proC <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],4],3][1]
        output <- rbind(output,c(complexinter$cell_a[i],complexinter$cell_b[i],proC,proD))}
    }else{
      proA <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],2],3][1]
      proB <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],3],3][1]
      proD <- complexinter$gene_a[i]
      output <- rbind(output,c(complexinter$cell_a[i],complexinter$cell_b[i],proD,proA),c(complexinter$cell_a[i],complexinter$cell_b[i],proD,proB))
      if (!complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],4] == ""){
        proC <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],4],3][1]
        output <- rbind(output,c(complexinter$cell_a[i],complexinter$cell_b[i],proD,proC))}
    }
  } 
}
list(lrpairs=unique(output), speed=speed, pairs=length(rownames(unique(output))))
}