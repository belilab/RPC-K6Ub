#load all necessary packages ----
library(DT)
library(plotly)
library(htmlwidgets)
library(openxlsx)
library(ggpubr)
library(ggrepel)
library(limma)
library(ggvenn)
library(ggrepel)
library(GGally)
library(biomaRt)
library(matrixTests)
library(ggnetwork)
library(network)
library(igraph)
library(clusterProfiler)
# define organism against the search/enrichment will be performed
org <- "org.Hs.eg.db"
library(org, character.only = TRUE)
library(RColorBrewer)
library(ReactomePA)
library(enrichplot)
library(ggraph)
library(tidygraph)
library(tidyverse)


### functions ----
`%!in%` = Negate(`%in%`)
#calculate significance with limma; define limma_function
limma_function = function(ratios, max_q = 0.05, min_count = 2, prefix = ''){
  result = data.frame(
    count = ncol(ratios) - rowSums(is.na(ratios)),
    mean = rowMeans(ratios, na.rm = T)
  )
  #fit = statistical testing using eBayes and lmFit to generate p and q (=FDR) values
  fit = as.data.frame(eBayes(lmFit(ratios[result$count >= min_count, ])))
  
  #adds the according p-value if the count value is more or equal the preset min_count
  result[result$count >= min_count, 'pvalue'] = fit$p.value
  #adds the accoridng q-value if the count value is more or equal the preset min_count; calculates the adjusted p-value with p.adjust using the fdr method
  #(here fdr is the same as BH = Benjamini & Hochberg correction; others can be chosen)
  result[result$count >= min_count, 'qvalue'] = p.adjust(fit$p.value, method = 'fdr')
  #up or douwn-regulation of a protein are computed giving a TRUE or FALSE statement
  results = mutate(
    result,
    up = count >= min_count & mean > 0 & qvalue <= max_q,
    down = count >= min_count & mean < 0 & qvalue <= max_q
  )
  names(results) = paste0(prefix, names(results))
  return(results)
}