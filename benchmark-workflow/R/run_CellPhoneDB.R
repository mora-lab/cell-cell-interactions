run_CellPhoneDB_output = function(name){
dir <- paste("Python/cpdb/",name,"/out/significant_means.txt", sep = "")
posi <- read.table(dir, header = T, sep = "\t")
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
for (i in 1:nrow(complexinter)) {
  if (complexinter[i,3] == "" & complexinter[i,4] == "") {
    proA <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],2],3]
    proB <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],3],3]
    proD <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],2],3]
    proE <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],3],3]
    output <- rbind(output,c(complexinter$cell_a[i],complexinter$cell_b[i],proA,proD),c(complexinter$cell_a[i],complexinter$cell_b[i],proA,proE),c(complexinter$cell_a[i],complexinter$cell_b[i],proB,proD),c(complexinter$cell_a[i],complexinter$cell_b[i],proB,proE))
    if (!complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],4] == ""){
      proC <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],4],3]
      output <- rbind(output,c(complexinter$cell_a[i],complexinter$cell_b[i],proC,proD),c(complexinter$cell_a[i],complexinter$cell_b[i],proC,proE))}
    if (!complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],4] == ""){
      proF <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],4],3]
      output <- rbind(output,c(complexinter$cell_a[i],complexinter$cell_b[i],proA,proF),c(complexinter$cell_a[i],complexinter$cell_b[i],proB,proF))}
    if (!complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],4] == "" & !complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],4] == ""){
      output <- rbind(output,c(complexinter$cell_a[i],complexinter$cell_b[i],proC,proF))}
  }
  else{
    if (complexinter[i,3] == "") {
      proA <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],2],3]
      proB <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],3],3]
      proD <- complexinter$gene_b[i]
      output <- rbind(output,c(complexinter$cell_a[i],complexinter$cell_b[i],proA,proD),c(complexinter$cell_a[i],complexinter$cell_b[i],proB,proD))
      if (!complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],4] == ""){
        proC <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_a[i],split = ":")[[1]][2],4],3]
        output <- rbind(output,c(complexinter$cell_a[i],complexinter$cell_b[i],proC,proD))}
    }
    else{
      proA <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],2],3]
      proB <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],3],3]
      proD <- complexinter$gene_a[i]
      output <- rbind(output,c(complexinter$cell_a[i],complexinter$cell_b[i],proD,proA),c(complexinter$cell_a[i],complexinter$cell_b[i],proD,proB))
      if (!complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],4] == ""){
        proC <- gene[gene$uniprot == complex[complex$complex_name == strsplit(complexinter$partner_b[i],split = ":")[[1]][2],4],3]
        output <- rbind(output,c(complexinter$cell_a[i],complexinter$cell_b[i],proD,proC))}
    }
  }  
}
print(paste("Total predicted L-R pairs:",nrow(output)))
output
}