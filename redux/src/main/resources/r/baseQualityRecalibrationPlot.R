library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

args <- commandArgs(trailing=T)
input <- args[1]
output <- gsub("tsv","png", input)

bqrDataRaw <- read.table(file = input, sep = "\t", header = T, stringsAsFactors = F)

COLORS6 = c("#2EBAED", "#000000", "#DE1C14", "#D4D2D2", "#ADCC54", "#F0D0CE")

standard_mutation <- function(types){
  types = gsub("G>T", "C>A", types)
  types = gsub("G>C", "C>G", types)
  types = gsub("G>A", "C>T", types)
  types = gsub("A>T", "T>A", types)
  types = gsub("A>G", "T>C", types)
  types = gsub("A>C", "T>G", types)
  return(types)
}

reverse_complement <- function(triContext) {
  complement = chartr("ATGC", "TACG", triContext)
  reverse = paste0(substring(complement, 3, 3), substring(complement, 2, 2), substring(complement, 1, 1))
  return (reverse)
}

prepare_data <- function(df, readType) {
  
  bqrData <- 
    df %>%
    filter(ref != alt) %>%
    mutate(
      changeInQual = recalibratedQual - originalQual,
      type = paste0(ref,">",alt), 
      substitution = factor(standard_mutation(type)),
      context = ifelse(type == substitution, trinucleotideContext, reverse_complement(trinucleotideContext)),
      context = factor(paste0(substring(context, 1, 1), ".", substring(context, 3, 3)))
    ) %>% 
    group_by(originalQual, substitution, context) %>%
    summarise(changeInQual = mean(changeInQual)) %>%
    ungroup() %>%
    spread(originalQual, changeInQual, fill = 0) %>%
    gather(originalQual, changeInQual, -1, -2) %>%
    select(substitution, context, originalQual, changeInQual) %>%
    mutate(originalQual = as.numeric(originalQual))
  
  uniqueBaseQuals = unique(bqrData$originalQual)
  if (length(uniqueBaseQuals) > 7) {
    sampleBy = ceiling(length(uniqueBaseQuals)/7)
    uniqueBaseQuals <- uniqueBaseQuals[seq(1, length(uniqueBaseQuals), sampleBy)]
    bqrData = bqrData %>% filter(originalQual %in% uniqueBaseQuals)
  }
  
  return(bqrData)
}

plot_data <- function(df) {
  ggplot(data = df, aes(x = context, y = changeInQual, fill = substitution, width = 0.6)) + 
    geom_bar(stat = "identity", colour = "black", size = 0.2) + 
    scale_fill_manual(values = COLORS6) + 
    facet_grid(originalQual ~ substitution) + 
    ylab("Base Quality Adjustment") + xlab("Context") + 
    guides(fill = FALSE) + theme_bw() +
    theme(axis.title.y = element_text(size = 12, vjust = 1), 
          axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
          axis.text.x = element_text(size = 5, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = 9),
          strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank())
}

bqrPlotsByReadType <- lapply(unique(bqrDataRaw$readType), function(readType) {
  
  bqrData <- prepare_data( bqrDataRaw[bqrDataRaw$readType==readType,] )
  bqrPlot <- plot_data(bqrData) + ggtitle(paste0("Read type: ", readType))
  
  return(bqrPlot)
})

bqrPlotsMerged <- patchwork::wrap_plots(bqrPlotsByReadType, ncol=1)

ggsave(
  filename = output, bqrPlotsMerged, 
  units = "in", height = 15, width = 12, scale = 1, dpi = 300
)
