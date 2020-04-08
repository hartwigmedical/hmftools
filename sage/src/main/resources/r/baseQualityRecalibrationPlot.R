library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(Biostrings)

COLORS6 = c("#2EBAED", "#000000", "#DE1C14", "#D4D2D2", "#ADCC54", "#F0D0CE")
SUBSTITUTIONS = c('C>A','C>G','C>T','T>A','T>C','T>G')
C_TRIPLETS = c("ACA", "ACC", "ACG", "ACT", "CCA", "CCC", "CCG", "CCT", "GCA", "GCC", "GCG", "GCT", "TCA", "TCC", "TCG", "TCT")
T_TRIPLETS = c( "ATA", "ATC", "ATG", "ATT", "CTA", "CTC", "CTG", "CTT", "GTA", "GTC", "GTG", "GTT", "TTA", "TTC", "TTG", "TTT")
CONTEXTS_96 = c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))

plot_96_profile_abs <- function (mut_matrix, colors)
{
  norm_mut_matrix = mut_matrix
  colors = COLORS6
  context = CONTEXTS_96
  substitution = rep(SUBSTITUTIONS, each = 16)
  substring(context, 2, 2) = "."
  df = data.frame(substitution = substitution, context = context)
  rownames(norm_mut_matrix) = NULL
  df2 = cbind(df, as.data.frame(norm_mut_matrix))
  df3 = melt(df2, id.vars = c("substitution", "context"))
  value = NULL

  plot = ggplot(data = df3, aes(x = context, y = value, fill = substitution, width = 0.6)) + 
    geom_bar(stat = "identity", colour = "black", size = 0.2) + 
    scale_fill_manual(values = colors) + 
    facet_grid(variable ~ substitution) + 
    ylab("Base Quality Adjustment") + xlab("Context") + 
    #coord_cartesian(ylim = c(0, ymax)) + scale_y_continuous(breaks = seq(0, ymax, 1)) + 
    guides(fill = FALSE) + theme_bw() + 
    theme(axis.title.y = element_text(size = 12, vjust = 1), 
          axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
          axis.text.x = element_text(size = 5, angle = 90, vjust = 0.4), 
          strip.text.x = element_text(size = 9),
          strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank())

  return(plot)
}

standard_mutation<-function(types){
  types = gsub("G>T", "C>A", types)
  types = gsub("G>C", "C>G", types)
  types = gsub("G>A", "C>T", types)
  types = gsub("A>T", "T>A", types)
  types = gsub("A>G", "T>C", types)
  types = gsub("A>C", "T>G", types)
  return(types)
}

args <- commandArgs(trailing=T)
#input <- "/Users/jon/hmf/analysis/COLO829T/sage/COLO829T.sage.bqr.tsv"
input <- args[1]
output <- gsub("tsv","png", input)

data = read.table(file = input, sep = "\t", header = T, stringsAsFactors = F) %>%
  filter(ref != alt) %>%
  mutate(
    label = originalQual,
    changeInQual = recalibratedQual - originalQual,
    type = paste0(ref,">",alt), 
    standardType = standard_mutation(type),
    standardContext = ifelse(type == standardType, trinucleotideContext, reverse(chartr("ATGC", "TACG", trinucleotideContext)))) %>%
  group_by(originalQual, standardType, standardContext) %>%
  summarise(changeInQual = mean(changeInQual)) %>%
  ungroup() %>%
  spread(originalQual, changeInQual) %>%
  arrange(
    standardType, 
    substring(standardContext, 1, 1), 
    substring(standardContext, 3, 3))

bqrPlot = plot_96_profile_abs(data.matrix(data %>% select(-1,-2)))
ggsave(filename = output, bqrPlot, units = "in", height = (dim(data)[2] - 1), width = 12, scale = 1)
