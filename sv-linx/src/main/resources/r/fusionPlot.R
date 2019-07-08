library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
theme_set(theme_bw())

# Parse the arguments
args <- commandArgs(trailing=T)
clusterProteinDomainPath <- args[1]
clusterFusedExonPath <- args[2]
circosPicturePath   <- args[3]

plot_fusion <- function(fusedExons, fusedProteinDomains) {
  
  fusion = fusedExons %>% group_by(fusion) %>% summarise(start = min(geneStart), end = max(geneEnd))
  
  fusedGenes = fusedExons %>% group_by(gene, upGene, color) %>% mutate(start = ifelse(rank == 1, start, geneStart)) %>%
    summarise(start = max(start, geneStart), end = max(geneEnd))
  
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

    geom_text(data = fusedGenes %>% filter(upGene), mapping = aes(label = gene, x = start, y = 1.1), hjust = 0, vjust = 0, size = 6 * 25.4 / 72) +
    geom_text(data = fusedGenes %>% filter(!upGene), mapping = aes(label = gene, x = end, y = 1.1), hjust = 1,  vjust = 0, size = 6 * 25.4 / 72) +
    scale_x_continuous(name = "", breaks = fusedExons$start, labels = fusedExons$rank) +
    scale_fill_manual(name =  "", values = proteinDomainColors) +
    scale_alpha_manual(values = fusedExonsAlpha) +
    coord_cartesian(ylim = c(0, 2)) +
    theme(axis.text.x = element_text(size = 5), axis.title = element_text(size = 5), legend.text = element_text(size = 5)) +
    theme(axis.ticks = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) +
    theme(panel.background = element_blank(), panel.border =  element_blank(), panel.grid = element_blank(), panel.spacing = unit(3, "pt")) +
    theme(plot.margin = margin(t = 0, b = 0, l = 3, r = 3, unit = "pt"), legend.box.margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "pt")) +
    theme(legend.position = c(0.5, 0.7), legend.key.size = unit(6, "pt"), legend.background=element_blank(), legend.key=element_blank(), legend.direction = "horizontal") +
    guides(fill = guide_legend(nrow = 1))
  
  if (nrow(fusedProteinDomains) > 0) {
    p1 = p1 + geom_rect(data = fusedProteinDomains, mapping = aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.5, fill = name), position = "identity", stat = "identity", alpha = 0.8)
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
for (selectedFusion in unique(clusterFusedExons$fusion)) {
  cat ("Processing", selectedFusion, "\n")

  fusedExons = clusterFusedExons %>%
    filter(fusion == selectedFusion) %>%
    mutate(
      upGene = ifelse(startsWith(fusion, gene) & rank == row_number(), T, F),
      color = ifelse(upGene, singleBlue, singleRed))
  
  fusedProteinDomains = clusterProteinDomains %>% 
    filter(chromosome == selectedFusion)
  
  fusionPlotList[[selectedFusion]] <- plot_fusion(fusedExons, fusedProteinDomains)
}

fusionPlotListCount = length(fusionPlotList)
pFusions = plot_grid(plotlist = fusionPlotList, ncol = 1)

pCircos <- ggdraw() + draw_image(circosPicturePath)

pWidth = 183
pHeight = pWidth * (1 + fusionPlotListCount * 0.07)

pCombined = plot_grid(pCircos, pFusions, ncol = 1, rel_heights = c(pWidth / pHeight, (pHeight - pWidth) / pWidth))
ggsave(circosPicturePath, pCombined, width = pWidth, height = pHeight, units = "mm", dpi = 300)



