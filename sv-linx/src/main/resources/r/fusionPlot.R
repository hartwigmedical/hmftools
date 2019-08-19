library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(magick)
theme_set(theme_bw())

clusterProteinDomainPath = "/Users/jon/hmf/analysis/fusions/data/CPCT02010990T.cluster54.COMPLEX.sv66.protein_domains.tsv"
clusterFusedExonPath = "/Users/jon/hmf/analysis/fusions/data/CPCT02010990T.cluster54.COMPLEX.sv66.fusions.tsv"
circosPicturePath = "/Users/jon/hmf/analysis/fusions/plot/CPCT02010990T.cluster54.COMPLEX.sv66.png"
fusionHeightPerRow = 250
fusionLegendHeightPerRow = 30
fontSize = 35

# Parse the arguments
args <- commandArgs(trailing=T)
clusterProteinDomainPath <- args[1]
clusterFusedExonPath <- args[2]
circosPicturePath   <- args[3]
fontSize <- as.numeric(args[4])

plot_fusion <- function(fusedExons, fusedProteinDomains, showLegend) {
  
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
  
  fusedExonsAlpha = setNames(c(1, fadedAlpha), c("false","true"))
  proteinDomainColors = fusedProteinDomains %>% select(name, color) %>% distinct()
  proteinDomainColors = setNames(proteinDomainColors$color, proteinDomainColors$name)
  
  p1 = ggplot() +
    geom_rect(data = fusion, mapping = aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.5), position = "identity", stat = "identity", fill = "#f5f5f5", color = NA) +
    geom_rect(data = fusedGenes, mapping = aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.5), position = "identity", stat = "identity", fill = fusedGenes$color, color = NA) +
    geom_rect(data = fusedExons, mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = 1), position = "identity", stat = "identity", fill = fusedExons$color, show.legend = F, color = NA) +

    geom_text(data = fusedGenes %>% filter(upGene), mapping = aes(label = label, x = start, y = 1.15), hjust = 0, vjust = 0, size = fontSize * 25.4 / 300) +
    geom_text(data = fusedGenes %>% filter(!upGene), mapping = aes(label = label, x = end, y = 1.15), hjust = 1,  vjust = 0, size = fontSize * 25.4 / 300) +
    scale_x_continuous(name = "", breaks = fusedExons$start, labels = fusedExons$rank) +
    scale_fill_manual(name =  "", values = proteinDomainColors) +
    scale_alpha_manual(values = fusedExonsAlpha) +
    coord_cartesian(ylim = c(0, 2)) +
    theme(axis.text.x = element_text(size = fontSize * 72 / 300), axis.title = element_text(size = 5), legend.text = element_text(size = 5)) +
    theme(axis.ticks = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) +
    theme(panel.background = element_blank(), panel.border =  element_blank(), panel.grid = element_blank(), panel.spacing = unit(3, "pt")) +
    theme(plot.margin = margin(t = 0, b = 0, l = 3, r = 3, unit = "pt"), legend.box.margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "pt"))
  
  if (nrow(fusedProteinDomains) > 0) {
    p1 = p1 + geom_rect(data = fusedProteinDomains, mapping = aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.5, fill = name), position = "identity", stat = "identity", alpha = 0.8)
  }
  
  # Want the faded area and segment break to appear on top of the protein domains
  p1 = p1 + 
    geom_rect(data = fadedArea, mapping = aes(xmin = start, xmax = end, ymin = 0.0, ymax = 1), position = "identity", stat = "identity", fill = "white", alpha = 1 - fadedAlpha, color = NA) + 
    geom_segment(data = fusedGenes %>% filter(upGene), mapping = aes(x = end, y = -0.1, xend = end, yend = 1.1))

  if (showLegend) {
    p1 = p1 +
      theme(legend.position = c(0.5, 0.8), 
            legend.key.size = unit(6, "pt"), 
            legend.text = element_text(size = fontSize * 72 / 300),
            legend.background=element_blank(), 
            legend.key=element_blank(), 
            legend.direction = "horizontal", 
            legend.spacing.x = unit(1, "mm")) +
      guides(fill = guide_legend(nrow = 1))
  } else {
    p1 = p1 + theme(legend.position = "none")
  }
  
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
  
  fusedProteinDomains = clusterProteinDomains %>% 
    filter(fusion == selectedFusion)
  
  pFusion = plot_fusion(fusedExons, fusedProteinDomains, T)
  pLegend = cowplot::get_legend(pFusion)
  pFusion = pFusion + theme(legend.position = "none")
  
  fusionPlotList[[selectedFusion]] <- pFusion
}

fusionPlotListCount = length(fusionPlotList)
fusionHeight = fusionHeightPerRow * fusionPlotListCount
fusionLegendHeight = fusionLegendHeightPerRow

pFusions = plot_grid(plotlist = fusionPlotList, ncol = 1)

imgCircos = magick::image_read(circosPicturePath)
pCircos <- ggdraw() + draw_image(imgCircos)

circosWidth = image_info(imgCircos)$width
circosHeight = image_info(imgCircos)$height

pCombined = plot_grid(pCircos, pFusions, pLegend, ncol = 1, rel_heights = c(circosHeight, fusionHeight, fusionLegendHeight))
ggsave(circosPicturePath, pCombined, width = circosWidth / 300, height = (circosHeight + fusionLegendHeight + fusionHeight) / 300, units = "in", dpi = 300, pointsize = 300)



