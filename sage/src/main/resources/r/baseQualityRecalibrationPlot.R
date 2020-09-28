#The MIT License (MIT)
#
#Copyright (c) 2016 Cuppen Research
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# Please note that the plotting logic was modified from MutationPatterns::plot_96_profile available from https://github.com/UMCUGenetics/MutationalPatterns/blob/master/R/plot_96_profile.R
#

library(ggplot2)
library(dplyr)
library(tidyr)

COLORS6 = c("#2EBAED", "#000000", "#DE1C14", "#D4D2D2", "#ADCC54", "#F0D0CE")

standard_mutation<-function(types){
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

prepare_data <- function(input) {
  read.table(file = input, sep = "\t", header = T, stringsAsFactors = F) %>%
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
}

plot_data <- function(df) {
  ggplot(data = df, aes(x = context, y = changeInQual, fill = substitution, width = 0.6)) + 
    geom_bar(stat = "identity", colour = "black", size = 0.2) + 
    scale_fill_manual(values = COLORS6) + 
    facet_grid(originalQual ~ substitution) + 
    ylab("Base Quality Adjustment") + xlab("Context") + 
    #coord_cartesian(ylim = c(0, ymax)) + scale_y_continuous(breaks = seq(0, ymax, 1)) + 
    guides(fill = FALSE) + theme_bw() + 
    theme(axis.title.y = element_text(size = 12, vjust = 1), 
          axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
          axis.text.x = element_text(size = 5, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = 9),
          strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank())
}

args <- commandArgs(trailing=T)
#input <- "/Users/jon/hmf/analysis/bqr/COLO829Tv001.sage.bqr.tsv"
input <- args[1]
output <- gsub("tsv","png", input)

bqrData = prepare_data(input)

uniqueBaseQuals = unique(bqrData$originalQual)
if (length(uniqueBaseQuals) > 7) {
  sampleBy = ceiling(length(uniqueBaseQuals)/7)
  uniqueBaseQuals <- uniqueBaseQuals[seq(1, length(uniqueBaseQuals), sampleBy)]
  bqrData = bqrData %>% filter(originalQual %in% uniqueBaseQuals)
}

bqrPlot = plot_data(bqrData)
ggsave(filename = output, bqrPlot, units = "in", height = length(unique(bqrData$originalQual)), width = 12, scale = 1, dpi = 300)
