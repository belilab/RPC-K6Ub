### GO-Enrichment using clusterProfiler 4.0 ----
# load df of Treamtnet1 vs UT with statistical analysis
df <- read.delim("limma_df_Treatment1_UT.txt", stringsAsFactors = F)

# filter for signifcantly upregulated proteins
temp_df <- df
temp_df <- subset(temp_df, temp_df$count >= 2 & temp_df$mean >= log2(2))
temp_df <- temp_df$Gene.names

# load total df as background to compare to
background <- unique(df$Gene.names)
universe <- background

# biomaRt (ENSEMBL)
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
  # View(datasets)
ensembl105 <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)

# ENSEMBL annotation of upregulated proteins
annotation_df <- getBM(attributes = c("entrezgene_id", "description", "external_gene_name"),
                       filters = c("external_gene_name"),
                       values = temp_df,
                       mart = ensembl105)

# ENSEMBL annotation of all quantified proteins
annotation_universe <- getBM(attributes = c("entrezgene_id", "description", "external_gene_name"),
                                filters = c("external_gene_name"),
                                values = universe,
                                mart = ensembl105)


##########################
### GO term annotation ###
##########################

# join both datasets, use the first ENSEMBL name for double annotated genes
#for high confidence hits
annotated_df <- temp_df
temp <- as.data.frame(matrix(nrow = length(temp_df), ncol=ncol(annotation_df)))
colnames(temp) <- colnames(annotation_df)
annotated_df <- cbind(annotated_df, temp)
colnames(annotated_df) <- c("Gene.names", colnames(temp))
annotated_df <- subset(annotated_df, annotated_df$Gene.names %in% annotation_df$external_gene_name)

for (i in 1:nrow(annotated_df)){
  temp <- subset(annotation_df, annotation_df$external_gene_name == annotated_df$Gene.names[i])
  if (nrow(temp) > 1){
    annotated_df$entrezgene_id[i] <- temp$entrezgene_id[1]
    annotated_df$description[i] <- temp$description[1]
    annotated_df$external_gene_name[i] <- temp$external_gene_name[1]
  }else{
    annotated_df$entrezgene_id[i] <- temp$entrezgene_id[1]
    annotated_df$description[i] <- temp$description[1]
    annotated_df$external_gene_name[i] <- temp$external_gene_name[1]
  }
}

# for all proteins
annotated_df_uni <- universe
temp <- as.data.frame(matrix(nrow = length(universe), ncol=ncol(annotation_universe)))
colnames(temp) <- colnames(annotation_universe)
annotated_df_uni <- cbind(annotated_df_uni, temp)
colnames(annotated_df_uni) <- c("Gene.names", colnames(temp))
annotated_df_uni <- subset(annotated_df_uni, annotated_df_uni$Gene.names %in% annotation_universe$external_gene_name)

for (i in 1:nrow(annotated_df_uni)){
  temp <- subset(annotation_universe, annotation_universe$external_gene_name == annotated_df_uni$Gene.names[i])
  if (nrow(temp) > 1){
    annotated_df_uni$entrezgene_id[i] <- temp$entrezgene_id[1]
    annotated_df_uni$description[i] <- temp$description[1]
    annotated_df_uni$external_gene_name[i] <- temp$external_gene_name[1]
  }else{
    annotated_df_uni$entrezgene_id[i] <- temp$entrezgene_id[1]
    annotated_df_uni$description[i] <- temp$description[1]
    annotated_df_uni$external_gene_name[i] <- temp$external_gene_name[1]
  }
}

# transform entrez_gene_id column into a vector and remove NA values
#for high confidence hits
ent_gene <- annotated_df$entrezgene_id
ent_gene <- ent_gene[which(!is.na(ent_gene))]
ent_gene <- as.character(ent_gene)
#for universe
ent_gene_uni <- annotated_df_uni$entrezgene_id
ent_gene_uni <- ent_gene_uni[which(!is.na(ent_gene_uni))]
ent_gene_uni <- as.character(ent_gene_uni)

# clusterProfiler GO enrichment - Molecular Function
ego_MF <- enrichGO(gene = ent_gene, OrgDb = org, ont = "MF", universe = ent_gene_uni,
                   pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)
write.table(x = ego_MF@results, file = "ego_MF.txt", sep = "\t")

# clusterProfiler GO enrichment - Biological Processes
ego_BP <- enrichGO(gene = ent_gene, OrgDb = org, ont = "BP", universe = ent_gene_uni,
                   pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)
write.table(x = ego_BP@results, file = "ego_BP.txt", sep = "\t")

# clusterProfiler GO enrichment - Cellular Compartments
ego_CC <- enrichGO(gene = ent_gene, OrgDb = org, ont = "CC", universe = ent_gene_uni,
                   pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)
write.table(x = ego_CC@results, file = "ego_CC.txt", sep = "\t")


######################################################################
# reading in data to generate lollipop plots
# filter data for Benjamini-Hochberg adjusted p-value <= 0.05
MF <- read.delim("ego_MF.txt", stringsAsFactors = FALSE)
MF <- subset(MF, p.adjust <= 0.05)
MF$origin <- "MF"
BP <- read.delim("ego_BP.txt", stringsAsFactors = FALSE)
BP <- subset(BP, p.adjust <= 0.05)
BP$origin <- "BP"
CC <- read.delim("ego_CC.txt", stringsAsFactors = FALSE)
CC <- subset(CC, p.adjust <= 0.05)
CC$origin <- "CC"

# assemble data in one common df
GO_data <- rbind(MF, BP, CC)

# calculate fold enrichment score
for (i in 1: nrow(GO_data)){
  GO_data$total[i] <- as.numeric(str_split(GO_data$GeneRatio[i], pattern = "/")[[1]][2])
  GO_data$total_bg[i] <- as.numeric(str_split(GO_data$BgRatio[i], pattern = "/")[[1]][2])
  GO_data$count_bg[i] <- as.numeric(str_split(GO_data$BgRatio[i], pattern = "/")[[1]][1])
  
  GO_data$fold_enrichment[i] <- (GO_data$Count[i]/GO_data$total[i])/(GO_data$count_bg[i]/GO_data$total_bg[i])
}

# order df in decreasing order according to fold enrichment
GO_data <- GO_data[order(GO_data$fold_enrichment, decreasing = TRUE),]

#reduce to top 25 terms from GO_data
GO_data <- GO_data[1:25, ]

# Plotting
levels_MF_ordered <- subset(GO_data, GO_data$origin == "MF")$Description
levels_BP_ordered <- subset(GO_data, GO_data$origin == "BP")$Description
levels_CC_ordered <- subset(GO_data, GO_data$origin == "CC")$Description
test_levels <- c(levels_BP_ordered, levels_MF_ordered, levels_CC_ordered)

GO_data <- within(GO_data,
                  Description <- factor(Description,
                                        levels = rev(test_levels)))

## plot
ggplot(data = GO_data, aes(x = Description, y = fold_enrichment)) +
  scale_x_discrete(expand = c(0.05,0.05),
                   breaks = waiver())+
  scale_y_continuous(expand = c(0,0.05),
                     breaks = waiver())+
  geom_segment(aes(x=Description,
                   xend=Description,
                   y=0,
                   yend=fold_enrichment),
               color="grey",
               linetype = "dashed") +
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_colour_gradient(low = "#CC1A1F", high = "#1472B9") +
  coord_flip() +
  theme(axis.line = element_line(colour = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.text = element_text(size = 12, face = "bold", colour = "black"),
        axis.title = element_text(size = 12, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "right") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("Fold enrichment") +
  xlab("top25 GO-terms")