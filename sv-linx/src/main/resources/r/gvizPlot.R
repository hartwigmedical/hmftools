library(Gviz)
library(dplyr)
library(grid)
library(ggplot2)
library(cowplot)

#chromosomeHeight = 150
#circosFontsize = 43
#circosPicturePath = "~/hmf/analysis/fusions/plot/JONTEST.png"
#bandsPath = "/Users/jon/hmf/analysis/fusions/data/cytoBand.txt"
#chromosomePath = "/Users/jon/hmf/analysis/fusions/data/COLO829T.cluster27.COMPLEX.sv8.chromosome.circos"

args <- commandArgs(trailing=T)
chromosomePath <- args[1]
bandsPath <- args[2]
circosPicturePath <- args[3]
circosFontsize <- as.numeric(args[4])
chromosomeHeight <- as.numeric(args[5])
  
doPlot <- function(contig) {
  start = chromosomeRanges[chromosomeRanges$chromosome == contig, "start"]
  end = chromosomeRanges[chromosomeRanges$chromosome == contig, "end"]
  chr = paste0("chr", contig)
  label = chromosomeRanges[chromosomeRanges$chromosome == contig, "label"]
  chrColor = chromosomeRanges[chromosomeRanges$chromosome == contig, "chrColor"]
  
  itrack <- IdeogramTrack(genome="X",chromosome=chr, name=label,bands=bands, ucscChromosomeNames=FALSE, 
                          fontsize = circosFontsize, fontcolor = "black", outline = T, min.height = 1, min.width = 1, fill = chrColor, col = chrColor, lwd = 10)
  plotTracks(itrack, from = start, to = end, add=TRUE, panel.only = T)
}

pCircos <- ggdraw() + draw_image(circosPicturePath)
bands <- read.table(bandsPath, h = T)
chromosomeRanges <- read.table(chromosomePath, h = F, sep = "\t", comment = "!", stringsAsFactors = F)
names(chromosomeRanges) <- c("chromosome", "start", "end", "chrColor")
chromosomeRanges = chromosomeRanges %>% mutate(label = paste0("chr ", chromosome))

chromosomeLengths = bands %>% 
  mutate(chromosome = gsub("chr", "", chrom)) %>% 
  group_by(chromosome) %>% 
  summarise(length = max(chromEnd)) %>% 
  filter(chromosome %in% chromosomeRanges$chromosome) %>%
  ungroup() %>%
  mutate(relLength = length / max(length)) %>%
  arrange(-relLength) %>%
  mutate(row = (row_number() - 1) %/% 4 + 1) %>%
  group_by(row) %>%
  mutate(rowLength = sum(relLength)) %>%
  ungroup() %>%
  mutate(width = relLength / max(rowLength)) %>%
  group_by(row) %>%
  mutate(
    priorWidth = lag(width),
    priorWidth = ifelse(is.na(priorWidth), 0, priorWidth),
    x = cumsum(priorWidth)) %>%
  ungroup() %>%
  mutate(height = 1 / max(row), y = (max(row) - row) * height) %>%
  select(chromosome, row, length, x, y, width, height)
chromosomeLengths = data.frame(chromosomeLengths)  

png(file = circosPicturePath, width = 3000, height = chromosomeHeight, units = "px")
for (i in 1:nrow(chromosomeLengths)) {
  x = chromosomeLengths[i, "x"]
  y = chromosomeLengths[i, "y"]
  width = chromosomeLengths[i, "width"]
  height = chromosomeLengths[i, "height"]
  contig = chromosomeLengths[i, "chromosome"]
   
  pushViewport(viewport(just = c("left","bottom"), x = x + 0.01, y = y, width = width - 0.02, height = height))
  doPlot(contig)
  popViewport(1)
}
dev.off()


pChr <- ggdraw() + draw_image(circosPicturePath)
pCombined = plot_grid(pChr, pCircos, ncol = 1, rel_heights = c(chromosomeHeight, 3000))
ggsave(circosPicturePath, pCombined, width = 3000/300.0, height = (3000 + chromosomeHeight)/300.0, units = "in", dpi = 300)
