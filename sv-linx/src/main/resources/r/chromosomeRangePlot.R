library(Gviz)
library(dplyr)
library(grid)
library(ggplot2)
library(cowplot)
library(magick)

#chromosomeHeightPerRow = 150
#chromosomeFontsize = 70
#chromosomeColumns = 4
#circosPicturePath = "~/hmf/analysis/fusions/plot/COLO829T.cluster27.COMPLEX.sv8.png"
#bandsPath = "/Users/jon/hmf/analysis/fusions/data/COLO829T.cluster27.COMPLEX.sv8.cytoBand.txt"
#chromosomeRangePath = "/Users/jon/hmf/analysis/fusions/data/COLO829T.cluster27.COMPLEX.sv8.chromosome.circos"

args <- commandArgs(trailing=T)
chromosomeRangePath <- args[1]
bandsPath <- args[2]
circosPicturePath <- args[3]
chromosomeFontsize <- as.numeric(args[4])
chromosomeHeightPerRow <- as.numeric(args[5])
chromosomeColumns <- as.numeric(args[6])

imgCircos = magick::image_read(circosPicturePath)
circosWidth = image_info(imgCircos)$width
  
doPlot <- function(contig) {
  start = chromosomeRanges[chromosomeRanges$chromosome == contig, "start"]
  end = chromosomeRanges[chromosomeRanges$chromosome == contig, "end"]
  chr = paste0("chr", contig)
  label = chromosomeRanges[chromosomeRanges$chromosome == contig, "label"]
  chrColor = chromosomeRanges[chromosomeRanges$chromosome == contig, "chrColor"]
  
  itrack <- IdeogramTrack(genome="X",chromosome=chr, name=label, bands=bands, ucscChromosomeNames=FALSE,
                          cex = 1, fontsize = chromosomeFontsize , fontcolor = "black", font = "helvetica",
                          fill = chrColor, col = chrColor, lwd = 10)
  plotTracks(itrack, from = start, to = end, add=TRUE, panel.only = F)
}

pCircos <- ggdraw() + draw_image(imgCircos)
bands <- read.table(bandsPath, h = T)
chromosomeRanges <- read.table(chromosomeRangePath, h = F, sep = "\t", comment = "!", stringsAsFactors = F)
names(chromosomeRanges) <- c("chromosome", "start", "end", "chrColor")
chromosomeRanges = chromosomeRanges %>% mutate(label = paste0("CHR ", chromosome))

chromosomeLengths = bands %>% 
  mutate(chromosome = gsub("chr", "", chrom)) %>% 
  group_by(chromosome) %>% 
  summarise(length = max(chromEnd)) %>% 
  filter(chromosome %in% chromosomeRanges$chromosome) %>%
  ungroup() %>%
  mutate(relLength = length / max(length)) %>%
  arrange(-relLength) %>%
  mutate(row = (row_number() - 1) %/% chromosomeColumns + 1) %>%
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

# Height is per row
chromosomeHeightPerRow = chromosomeHeightPerRow * max(chromosomeLengths$row)

#circosPicturePath = "~/hmf/analysis/fusions/plot/COLO829T.cluster27.COMPLEX.sv8.out.png"
png(file = circosPicturePath, width = circosWidth, height = chromosomeHeightPerRow, units = "px")
for (i in 1:nrow(chromosomeLengths)) {
  x = chromosomeLengths[i, "x"]
  y = chromosomeLengths[i, "y"]
  width = chromosomeLengths[i, "width"]
  height = chromosomeLengths[i, "height"]
  contig = chromosomeLengths[i, "chromosome"]
   
  pushViewport(viewport(just = c("left","bottom"), x = x, y = y, width = width, height = height))
  doPlot(contig)
  popViewport(1)
}
dev.off()


pChr <- ggdraw() + draw_image(circosPicturePath)
pCombined = plot_grid(pCircos, pChr, ncol = 1, rel_heights = c(circosWidth, chromosomeHeightPerRow - 0.1 * chromosomeHeightPerRow))
ggsave(circosPicturePath, pCombined, width = circosWidth/300.0, height = (circosWidth + chromosomeHeightPerRow)/300.0, units = "in", dpi = 300)
