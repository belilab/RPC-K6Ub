# Volcano Plot
# subsetting for diGly sites with at least 2 quantifications
pg_gg <- subset(pg_gg, pg_gg$count >= 2)

#dot color and labeling according to significance cutoffs (qvalue <= 0.05 & log2(FC) >= 2)
for (i in 1: nrow(pg_gg)){
  if (pg_gg$mean[i] >= 2 & pg_gg$qvalue[i] <= 0.05)
  {pg_gg$threshold[i] <- "#cb181d"
  pg_gg$label[i] <- 1}
  else if (pg_gg$mean[i] >= 1 & pg_gg$mean[i] <= 2 & pg_gg$qvalue[i] <= 0.05)
  {pg_gg$threshold[i] <- "#fc9272"
  pg_gg$label[i] <- 2}
  else if (pg_gg$mean[i] <= -2 & pg_gg$qvalue[i] <= 0.05)
  {pg_gg$threshold[i] <- "#0073C2FF"
  pg_gg$label[i] <- 3}
  else if (pg_gg$mean[i] <= -1 & pg_gg$mean[i] >= -2 & pg_gg$qvalue[i] <= 0.05)
  {pg_gg$threshold[i] <- "#9ecae1"
  pg_gg$label[i] <- 4}
  else {pg_gg$threshold[i] <- "#bdbdbd"
  pg_gg$label[i] <- 0}
}

#Plotting

#plot with ggplot
volcano <- ggplot(data=pg_gg, aes(x = mean, y = -log10(qvalue), 
                                  text = paste("Protein:", Gene.names, "\n",
                                               "GlyGly..K..Probabilities", GlyGly..K..Probabilities, "\n",
                                               "fold-change:", ifelse(mean > 0, round(2^abs(mean), digits = 3), -round(2^abs(mean), digits = 3)),"\n", 
                                               "-log10(q-value):", round(-log10(qvalue), digits = 3), "\n",
                                               "Description:", Protein.names, "\n",
                                               "K-Position:", Position, "\n",
                                               sep = " "))) +
  scale_x_continuous(name="log2(Treatment1 / UT)",
                     expand = c(0,0),
                     limits = c(min(pg_gg$mean)-0.5, max(pg_gg$mean)+0.5),
                     breaks = waiver(),
                     n.breaks = 5)+
  scale_y_continuous(name="-log10(q-value)",
                     expand = c(0,0),
                     limits = c(0, -log10(min(pg_gg$qvalue))+0.5),
                     breaks = waiver(),
                     n.breaks = 8)+
  geom_hline(yintercept=-log10(0.05),
             linetype="dashed",
             color = "darkgrey",
             size = 0.7,
             alpha = 0.7) +
  geom_vline(xintercept=2,
             linetype="dashed",
             color = "darkgrey",
             size = 0.7, 
             alpha = 0.7) +
  geom_vline(xintercept=1,
             linetype="dashed",
             color = "lightgrey",
             size = 0.7) +
  geom_vline(xintercept=-2,
             linetype="dashed",
             color = "darkgrey",
             size = 0.7,
             alpha = 0.7) +
  geom_vline(xintercept=-1,
             linetype="dashed",
             color = "lightgrey",
             size = 0.7) +
  geom_point(aes(alpha=1,
                 colour=threshold),
             size = 3) +
  geom_point(data =filter(pg_gg, label == 1),
             shape = 1,
             size = 3,
             colour = "black")+
  geom_point(data =filter(pg_gg, label == 3),
             shape = 1,
             size = 3,
             colour = "black")+
  geom_label_repel(data=filter(pg_gg, label == 1),
                   aes(label=paste0(Gene.names, " (", Amino.acid, Position, ")")),
                   na.rm = TRUE,
                   max.overlaps = 500,
                   max.time = 2,
                   size = 4,
                   segment.colour = 'black',
                   segment.alpha = 0.5,
                   # hjust = -0.1,
                   segment.curvature = 0.1,
                   segment.inflect = TRUE,
                   segment.square = FALSE,
                   xlim = c(0, NA)) +
  scale_color_identity()+
  theme(axis.line = element_line(colour = "black", size = 1),
        axis.ticks = element_line(colour = "black", size = 1),
        axis.text = element_text(size = 14, face = "bold", colour = "black"),
        axis.title = element_text(size = 14, face = "bold", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour = "black", size = 1),
        legend.position="none",
        panel.background = element_blank(),) +
  theme(plot.title = element_text(hjust = 0.5))

#creating output pdf file
pdf("working_directory/volcanoplot_Treatment1vsUT.pdf", height = 13, width = 14)
volcano
dev.off()