### PREPROCESSING ----
#set working directory
setwd(dir ="working_directory")

#read data of replicate 1-3
pg_gg <- read.delim("GlyGly (K)Sites.txt", stringsAsFactors = FALSE)
raw <- pg_gg

#set working directory to analysis folder
setwd(dir ="working_directory/analysis")

#reduce regular Gene.names strings
pg_gg$Gene.names <- sub(";.*$", "", pg_gg$Gene.names)

#define subset in proteinGroups; filtering for reverse positive, contaminant, and diGly-modification localization probability >= 0.95
no_reverse <-pg_gg [,"Reverse"] != "+"
no_contaminant <- pg_gg [,"Potential.contaminant"] != "+"
low_localization_prob <- pg_gg [, "Localization.prob"] >= 0.95

#filter to subset
pg_gg <- subset(pg_gg, no_reverse & no_contaminant & low_localization_prob)

rm(no_contaminant,no_reverse, low_localization_prob)

#log2 transformation and define ratios
# # log2 transformation of MaxQuant-normalized ratios and Label switch (SILAC)
pg_gg$log2.Treatment1.UT.1 <- log2(pg_gg$Ratio.H.L.normalized.1)
pg_gg$log2.Treatment1.UT.2 <- log2(1/pg_gg$Ratio.M.L.normalized.2)
pg_gg$log2.Treatment1.UT.3 <- log2(pg_gg$Ratio.H.L.normalized.3)


pg_gg$log2.Treatment2.UT.1 <- log2(pg_gg$Ratio.M.L.normalized.1)
pg_gg$log2.Treatment2.UT.2 <- log2(pg_gg$Ratio.H.M.normalized.2)
pg_gg$log2.Treatment2.UT.3 <- log2(pg_gg$Ratio.M.L.normalized.3)


pg_gg$log2.Treatment1.Treatment2.1 <- log2(pg_gg$Ratio.H.M.normalized.1)
pg_gg$log2.Treatment1.Treatment2.2 <- log2(1/pg_gg$Ratio.H.L.normalized.2)
pg_gg$log2.Treatment1.Treatment2.3 <- log2(pg_gg$Ratio.H.M.normalized.3)


#add missing gene names from FASTA header
missing_gene_names <- which(pg_gg$Gene.names == "")
for (i in missing_gene_names){
  #adding gene.names
  pg_gg$Gene.names[i] <- sub(" PE=.*", "", pg_gg$Fasta.headers[i])
  pg_gg$Gene.names[i] <- sub(".*GN=","",pg_gg$Gene.names[i])
  #adding protein names/description
  pg_gg$Protein.names[i] <- sub(" OS=.*", "", pg_gg$Fasta.headers[i])
  pg_gg$Protein.names[i] <- sub(".*HUMAN ","",pg_gg$Protein.names[i])
  
  if (startsWith(x = pg_gg$Gene.names[i], prefix = ">")){
    pg_gg$Gene.names[i] <- pg_gg$Protein.names[i]
  }
}

# reducing dataframe
pg_gg_filtered <- pg_gg
pg_gg <- pg_gg_filtered[ ,c(5:6, 8:13, 27:29, 32:33, 265, 282:290)]
