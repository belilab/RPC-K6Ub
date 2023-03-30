## LIMMA for normalized Log2- transformed ratios

# run the limma_function defined in Packages_and_functions.R
#define experiment names
Exp_name <- c("log2.Treatment1.UT.1", "log2.Treatment1.UT.2", "log2.Treatment1.UT.3")
#calculate statistics for diGly-sites that were at least 2 times quantified via MS
hl_statistics = limma_function(pg_gg[ ,Exp_name], min_count = 2)
#add statistical analysis to working df
pg_gg = bind_cols(pg_gg, hl_statistics)
pg_gg$Signature <- paste0(pg_gg$Gene.names, "_", pg_gg$Amino.acid, pg_gg$Position)
#selecting relevant columns and creating output df
output <- pg_gg[ ,c(1:2, 30, 3:17, 24:29)]
#writing output file
write.table(x = output, file = "limma_Treatment1vsUT.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)