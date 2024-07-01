library(Gviz)
library(dplyr)
library(grid)
library(ggplot2)
library(cowplot)
library(magick)


args <- commandArgs(trailing=T)
chromosomeRangePath <- args[1]
bandsPath <- args[2]
circosPicturePath <- args[3]
chromosomeFontsize <- as.numeric(args[4])
chromosomeHeightPerRow <- as.numeric(args[5])
chromosomeMaxColumns <- as.numeric(args[6])

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
chromosomeRanges = chromosomeRanges %>% mutate(chromosome = gsub("chr", "", chromosome), label = paste0("CHR ", chromosome))

add_row_numbers <- function(columnsInFirstRow, chromosomeLengths, minRelativeLength = 0.1) {
  lengthOfFirstRow = sum(chromosomeLengths %>% filter(row_number() <= columnsInFirstRow) %>% pull(length))
  chromosomeLengths = chromosomeLengths %>% 
    mutate(
      relLength = pmax(minRelativeLength, length / lengthOfFirstRow), 
      totalLength = cumsum(relLength), 
      row = ceiling(totalLength))
  
  return (chromosomeLengths)
}

optimal_columns_in_first_row <- function(chromosomeMaxColumns, chromosomeMaxRows, minRelativeLength, chromosomeLengths) {
  for(columns in 1:chromosomeMaxColumns) {
    tmp = add_row_numbers(columns, chromosomeLengths, minRelativeLength)
    rows = max(tmp$row)
    if(rows <= chromosomeMaxRows) {
      return (columns)
    }
  }
  return (optimal_columns_in_first_row(chromosomeMaxColumns, rows, minRelativeLength, chromosomeLengths))  
}

chromosomeLengths = bands %>% 
  mutate(
    chromosome = gsub("chr", "", chrom), 
    chromosome = factor(chromosome, levels = c(1:22, "X", "Y"), ordered = T)) %>% 
  group_by(chromosome) %>% 
  summarise(length = max(chromEnd)) %>% 
  filter(chromosome %in% chromosomeRanges$chromosome) %>%
  ungroup() %>%
  arrange(chromosome) %>%
  mutate(relLength = length / max(length))

chromosomeMaxRows = ceiling(nrow(chromosomeRanges) / chromosomeMaxColumns)
optimalColumsInFirstRow = optimal_columns_in_first_row(chromosomeMaxColumns, chromosomeMaxRows, 0.1, chromosomeLengths)
chromosomeLengths = add_row_numbers(optimalColumsInFirstRow, chromosomeLengths, 0.1)

chromosomeLengths = chromosomeLengths %>%
  group_by(row) %>%
  mutate(
    column = row_number(), 
    rowLength = sum(relLength),
    spacing = (1 - rowLength) / 2.0,
    spacing = ifelse(column == 1, spacing, 0)) %>%
  ungroup() %>%
  mutate(width = relLength / max(rowLength)) %>%
  group_by(row) %>%
  mutate(
    priorWidth = lag(width),
    priorWidth = ifelse(is.na(priorWidth), 0, priorWidth),
    x = cumsum(priorWidth) + cumsum(spacing)) %>%
  ungroup() %>%
  mutate(height = 1 / max(row), y = (max(row) - row) * height) %>%
  select(chromosome, row, length, x, y, width, height)
chromosomeLengths = data.frame(chromosomeLengths)

# Height is per row
chromosomeHeight = chromosomeHeightPerRow * max(chromosomeLengths$row)

png(file = circosPicturePath, width = 0.85 * circosWidth, height = chromosomeHeight, units = "px")

for(i in 1:nrow(chromosomeLengths)) {
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
pChr = pChr + theme(plot.margin = margin(t = 0, b = 0, l = 3, r = 3, unit = "pt"))

pCombined = plot_grid(pCircos, pChr, ncol = 1, rel_heights = c(circosWidth, chromosomeHeight))
ggsave(circosPicturePath, pCombined, width = circosWidth/300.0, height = (circosWidth + chromosomeHeight)/300.0, units = "in", dpi = 300, bg="white")
