stringsAsFactors=FALSE

inter <- read.csv("interaction_input.csv")
gene <- read.csv("gene_input.csv")
complex <- read.csv("complex_input.csv")


singleinter <- inter[!(inter$protein_name_a == ""|inter$protein_name_b == ""),]
db <- NULL
for (i in 1:nrow(singleinter)) {
db$geneA[i] <- gene[gene$uniprot == singleinter$partner_a[i],3]
db$geneB[i] <- gene[gene$uniprot == singleinter$partner_b[i],3]
}
db <- as.data.frame(db)
db <- db[!(db$geneA %in% "" | db$geneB %in% ""),]


complexinter <- inter[inter$protein_name_a == ""|inter$protein_name_b == "",]

for (i in 1:nrow(complexinter)) {
  if (complexinter[i,4] == "" & complexinter[i,5] == "") {
    proA <- gene[gene$uniprot == complex[complex$complex_name == complexinter$partner_a[i],2],3]
    proB <- gene[gene$uniprot == complex[complex$complex_name == complexinter$partner_a[i],3],3]
    proD <- gene[gene$uniprot == complex[complex$complex_name == complexinter$partner_b[i],2],3]
    proE <- gene[gene$uniprot == complex[complex$complex_name == complexinter$partner_b[i],3],3]
    db <- rbind(db,c(proA,proD),c(proA,proE),c(proB,proD),c(proB,proE))
    if (!complex[complex$complex_name == complexinter$partner_a[i],4] == ""){
      proC <- gene[gene$uniprot == complex[complex$complex_name == complexinter$partner_a[i],4],3]
      db <- rbind(db,c(proC,proD),c(proC,proE))}
    if (!complex[complex$complex_name == complexinter$partner_b[i],4] == ""){
      proF <- gene[gene$uniprot == complex[complex$complex_name == complexinter$partner_b[i],4],3]
      db <- rbind(db,c(proA,proF),c(proB,proF))}
    if (!complex[complex$complex_name == complexinter$partner_a[i],4] == "" & !complex[complex$complex_name == complexinter$partner_b[i],4] == ""){
      db <- rbind(db,c(proC,proF))}
  }
  else{
  if (complexinter[i,4] == "") {
    proA <- gene[gene$uniprot == complex[complex$complex_name == complexinter$partner_a[i],2],3]
    proB <- gene[gene$uniprot == complex[complex$complex_name == complexinter$partner_a[i],3],3]
    proD <- strsplit(complexinter[i,5],split="_")[[1]][1]
    db <- rbind(db,c(proA,proD),c(proB,proD))
    if (!complex[complex$complex_name == complexinter$partner_a[i],4] == ""){
     proC <- gene[gene$uniprot == complex[complex$complex_name == complexinter$partner_a[i],4],3]
     db <- rbind(db,c(proC,proD))}
  }
  else{
    proA <- gene[gene$uniprot == complex[complex$complex_name == complexinter$partner_b[i],2],3]
    proB <- gene[gene$uniprot == complex[complex$complex_name == complexinter$partner_b[i],3],3]
    proD <- strsplit(complexinter[i,4],split="_")[[1]][1]
    db <- rbind(db,c(proD,proA),c(proD,proB))
    if (!complex[complex$complex_name == complexinter$partner_b[i],4] == ""){
      proC <- gene[gene$uniprot == complex[complex$complex_name == complexinter$partner_b[i],4],3]
      db <- rbind(db,c(proD,proC))}
  }
  }  
}

write.table(db,"CellPhoneDB database 2.0.0.txt", sep = "\t", quote = F, row.names = F)
