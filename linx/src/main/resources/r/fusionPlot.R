library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(magick)
theme_set(theme_bw())

# Parse the arguments
args <- commandArgs(trailing=T)
clusterProteinDomainPath <- args[1]
clusterFusedExonPath <- args[2]
circosPicturePath   <- args[3]
fontSize <- as.numeric(args[4])
fusionLegendRows <- as.numeric(args[5])
fusionLegendHeightPerRow <- as.numeric(args[6])
fusionHeightPerRow <- as.numeric(args[7])

plot_legend <- function(nonUTRProteinDomain) {
  
  proteinDomainColors = nonUTRProteinDomain %>% select(name, color) %>% distinct()
  proteinDomainColors = setNames(proteinDomainColors$color, proteinDomainColors$name)
  
  p1 <- ggplot(nonUTRProteinDomain ) + 
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = name)) +
    theme(legend.position = "bottom", legend.key.size = unit(6, "pt"), legend.text = element_text(size = fontSize * 72 / 300), legend.background=element_blank(), legend.key=element_blank(), legend.direction = "horizontal", legend.spacing.x = unit(1, "mm")) +
    guides(fill = guide_legend(nrow = fusionLegendRows)) + 
    scale_fill_manual(name =  "", values = proteinDomainColors, drop = F)
  return (p1)  
}

plot_fusion <- function(fusedExons, fusedProteinDomains) {
  
  gene = fusedExons %>% group_by(transcript, color) %>% summarise()
  utrRegions = fusedProteinDomains %>% filter(name == "UTR/Non-coding")
  otherRegions = fusedProteinDomains %>% filter(name != "UTR/Non-coding") %>% arrange(name)
  
  fusion = fusedExons %>% group_by(fusion) %>% summarise(start = min(geneStart), end = max(geneEnd))
  
  fusedGenes = fusedExons %>% 
    group_by(gene, transcript, upGene, color) %>% 
    mutate(start = ifelse(rank == 1, start, geneStart)) %>%
    summarise(start = max(start, geneStart), end = max(geneEnd)) %>%
    mutate(label = paste0(gene, " - ",transcript))
  
  fadedAlpha = 0.4
  fadedArea = fusedExons %>% 
    filter(skipped == "false") %>% 
    group_by(fusion, upGene) %>% 
    summarise(start = min(start), end = max(end)) %>%
    mutate(value = ifelse(upGene, end, start), upGene = ifelse(upGene, "start", "end")) %>% 
    select(-start, -end)   %>% spread(upGene, value)
  
  proteinDomainColors = fusedProteinDomains %>% filter(name != "UTR/Non-coding") %>% select(name, color) %>% distinct()
  proteinDomainColors = setNames(proteinDomainColors$color, proteinDomainColors$name)
  
  p1 = ggplot() +
    geom_rect(data = fusedGenes, mapping = aes(xmin = start, xmax = end, ymin = 0.45, ymax = 0.55), position = "identity", stat = "identity", fill = fusedGenes$color, color = NA) +
    geom_rect(data = fusedExons, mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = 1), position = "identity", stat = "identity", fill = fusedExons$color, show.legend = F, color = NA) +
    geom_text(data = fusedGenes %>% filter(upGene), mapping = aes(label = label, x = start, y = 1.15), hjust = 0, vjust = 0, size = fontSize * 25.4 / 300) +
    geom_text(data = fusedGenes %>% filter(!upGene), mapping = aes(label = label, x = end, y = 1.15), hjust = 1,  vjust = 0, size = fontSize * 25.4 / 300) +
    scale_x_continuous(name = "", breaks = fusedExons$start, labels = fusedExons$rank) +
    scale_fill_manual(name =  "", values = proteinDomainColors) +
    coord_cartesian(ylim = c(0, 2)) +
    theme(axis.text.x = element_text(size = fontSize * 72 / 300), axis.title = element_text(size = 5), legend.text = element_text(size = 5)) +
    theme(axis.ticks = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) +
    theme(panel.background = element_blank(), panel.border =  element_blank(), panel.grid = element_blank(), panel.spacing = unit(3, "pt")) +
    theme(plot.margin = margin(t = 0, b = 0, l = 3, r = 3, unit = "pt"), legend.box.margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "pt")) + 
    theme(legend.position = "none")

  if (nrow(utrRegions) > 0) {
    p1 = p1 + 
      geom_rect(data = utrRegions, mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = 1), position = "identity", stat = "identity", alpha = 1, fill = "white") +
      geom_rect(data = utrRegions, mapping = aes(xmin = start, xmax = end, ymin = 0.3, ymax = 0.7), position = "identity", stat = "identity", fill = utrRegions$color)
  }
  
  if (nrow(otherRegions) > 0) {
    p1 = p1 + geom_rect(data = otherRegions, mapping = aes(xmin = start, xmax = end, ymin = 0.3, ymax = 0.7, fill = name), position = "identity", stat = "identity", alpha = 0.8)
  }
  
  # Want the faded area and segment break to appear on top of the protein domains
  p1 = p1 + 
    geom_rect(data = fadedArea, mapping = aes(xmin = start, xmax = end, ymin = 0.0, ymax = 1), position = "identity", stat = "identity", fill = "white", alpha = 1 - fadedAlpha, color = NA) + 
    geom_segment(data = fusedGenes %>% filter(upGene), mapping = aes(x = end, y = -0.1, xend = end, yend = 1.1))


  return (p1)
}

singleBlue = "#6baed6"
singleRed = "#d6906b"

clusterProteinDomains = read.table(clusterProteinDomainPath, sep = '\t', header = T, comment.char = "$", stringsAsFactors = F) %>%
  mutate(name = gsub(" domain", "", name))
clusterFusedExons = read.table(clusterFusedExonPath, sep = '\t', header = T, stringsAsFactors = F)

fusionPlotList = list()
fusions = unique(clusterFusedExons$fusion)

for (i in c(1:length(fusions))) {
  selectedFusion = fusions[i]
  cat ("Processing", selectedFusion, "\n")

  fusedExons = clusterFusedExons %>%
    filter(fusion == selectedFusion) %>%
    mutate(
      upGene = ifelse(startsWith(fusion, gene) & rank == row_number(), T, F),
      color = ifelse(upGene, singleBlue, singleRed))
  
  fusedProteinDomains = clusterProteinDomains %>% filter(fusion == selectedFusion)
  
  pFusion = plot_fusion(fusedExons, fusedProteinDomains)
  fusionPlotList[[selectedFusion]] <- pFusion
}

nonUTRProteinDomain = clusterProteinDomains %>% filter(name != "UTR/Non-coding")
if (nrow(nonUTRProteinDomain) > 0) {
  pLegend = plot_legend(nonUTRProteinDomain)
  pLegend = cowplot::get_legend(pLegend)
} else {
  pLegend = NA
}

fusionPlotListCount = length(fusionPlotList)
fusionHeight = fusionHeightPerRow * fusionPlotListCount
fusionLegendHeight = fusionLegendHeightPerRow * fusionLegendRows

pFusions = plot_grid(plotlist = fusionPlotList, ncol = 1)

imgCircos = magick::image_read(circosPicturePath)
pCircos <- ggdraw() + draw_image(imgCircos)

circosWidth = image_info(imgCircos)$width
circosHeight = image_info(imgCircos)$height

if (is.na(pLegend)) {
  totalHeight = (circosHeight + fusionHeight) 
  pCombined = plot_grid(pCircos, pFusions, ncol = 1, rel_heights = c(circosHeight, fusionHeight))
} else {
  totalHeight = (circosHeight + fusionLegendHeight + fusionHeight) 
  pCombined = plot_grid(pCircos, pFusions, pLegend, ncol = 1, rel_heights = c(circosHeight, fusionHeight, fusionLegendHeight))
}
ggsave(circosPicturePath, pCombined, width = circosWidth / 300, height =  totalHeight/ 300, units = "in", dpi = 300, pointsize = 300)



