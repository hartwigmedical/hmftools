library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
theme_set(theme_bw())

singleBlue = "#6baed6"
singleRed = "#d94701"


proteinDomain = read.table("/Users/jon/hmf/analysis/CPCT02390005T/data/CPCT02390005T.cluster088.COMPLEX.sv18.protein_domains.tsv", sep = '\t', header = T, comment.char = "$", stringsAsFactors = F) %>%
  mutate(name = gsub("domain", "", name))
proteinDomainColors = proteinDomain %>% select(name, color) %>% distinct()
proteinDomainColors = setNames(proteinDomainColors$color, proteinDomainColors$name)
fusedExons = read.table("/Users/jon/hmf/analysis/CPCT02390005T/data/CPCT02390005T.cluster088.COMPLEX.sv18.fusions.tsv", sep = '\t', header = T, stringsAsFactors = F) %>% 
  mutate(upGene = ifelse(startsWith(fusion, gene), T, F), color = ifelse(upGene, "#6baed6", "#d6906b"))
fusedGene = fusedExons %>% select(gene, start = geneStart, end = geneEnd, upGene, color) %>% distinct()
fusedExonsAlpha = setNames(c(1, 0.4), c("false","true"))

ggplot() +
  geom_rect(data = fusedGene, mapping = aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.5), position = "identity", stat = "identity", fill = fusedGene$color) + 
  geom_rect(data = fusedExons, mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = 1, alpha = truncated), position = "identity", stat = "identity", fill = fusedExons$color, show.legend = F) +
  geom_rect(data = proteinDomain, mapping = aes(xmin = start, xmax = end, ymin = 0.0, ymax = 0.5, fill = name), position = "identity", stat = "identity", alpha = 0.8) +
  geom_text(data = fusedGene %>% filter(upGene), mapping = aes(label = gene, x = start, y = 1.1), hjust = 0, vjust = 0, size = 8 * 25.4 / 72) +
  geom_text(data = fusedGene %>% filter(!upGene), mapping = aes(label = gene, x = end, y = 1.1), hjust = 1,  vjust = 0, size = 8 * 25.4 / 72) +
  scale_x_continuous(name = "", breaks = fusedExons$start, labels = fusedExons$rank) +
  scale_fill_manual(name =  "", values = proteinDomainColors) +
  scale_alpha_manual(values = fusedExonsAlpha) +
  coord_cartesian(ylim = c(0, 1.8)) +
  theme(axis.ticks = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.text.x = element_text(size = 7), axis.title = element_text(size = 7)) +
  theme(panel.background = element_blank(), panel.border =  element_blank(), panel.grid = element_blank(), panel.spacing = unit(3, "pt")) +
  theme(legend.text = element_text(size = 7), legend.position = c(0.5, 0.8), legend.margin = margin(t = 0, b = 0, l = 3, r = 3, unit = "pt")) +
  theme(plot.margin = margin(t = 0, b = 0, l = 3, r = 3, unit = "pt"), legend.box.margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "pt")) +
  theme(legend.key.size = unit(10, "pt"), legend.background=element_blank(), legend.key=element_blank(), legend.direction = "horizontal")

  


ggplot() +

    geom_rect(data = geneDetails, mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5), position = "identity", stat = "identity", fill= "blue", alpha = 0.3) +
  geom_rect(data = geneExons, mapping = aes(xmin = start, xmax = end, ymin = 0.5, ymax = 0.7), position = "identity", stat = "identity", fill= "blue", alpha = 0.3) +
  geom_rect(data = geneExons, mapping = aes(xmin = start, xmax = end, ymin = -0.2, ymax = 0), position = "identity", stat = "identity", fill= "blue", alpha = 0.3) +
  geom_rect(data = geneProteinDomain, mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5), position = "identity", stat = "identity", fill= "red", alpha = 0.8) +
  geom_text(data = geneProteinDomain, mapping = aes(x = start, y = 0.5, label = name, hjust = 0, vjust = 1)) +
  scale_x_continuous(name = "Exon Rank", breaks = geneExons$start, labels = geneExons$rank) +
  theme(axis.ticks = element_blank(), panel.background = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) +
  ggtitle(geneName)



exons = read.table("/Users/jon/hmf/analysis/CPCT02390005T/data/CPCT02390005T.cluster088.COMPLEX.sv18.exon.circos") %>% 
  separate(V5, c("gene", "rank"),sep=",") %>%
  mutate(
    rank = substring(rank, 6),
    gene = substring(gene, 6)) %>%
  group_by(gene) %>%
  mutate(
    chromosome = substring(V1, 3),
    start = V2 - min(V2),
    end = V3 - min(V2)) 

proteinDomain = read.table("/Users/jon/hmf/analysis/CPCT02390005T/data/CPCT02390005T.cluster088.COMPLEX.sv18.protein_domain.circos") %>% 
  separate(V5, c("fill", "color", "name"),sep=",") %>%
  mutate(
    chromosome = substring(V1, 3),
    fill = substring(fill, 12),
    color = substring(color, 7),
    name = substring(name, 6)
  )


geneName = "KIF5B"

geneDetails = exons %>% 
  filter(gene == geneName) %>%
  group_by(gene, chromosome) %>%
  summarise(V2 = min(V2), V3 = max(V3), start = min(start), end = max(end))

geneChromosome = geneDetails$chromosome
geneMin = geneDetails$V2
geneMax = geneDetails$V3

geneExons = exons %>% filter(gene == geneName) %>% mutate(middle = (end - start) / 2)

geneProteinDomain = proteinDomain %>% filter(chromosome == geneChromosome, V2 >= geneMin, V3 <= geneMax) %>%
  mutate(
    start = V2 - geneMin,
    end = V3 - geneMin)

ggplot() +
  geom_rect(data = geneDetails, mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5), position = "identity", stat = "identity", fill= "blue", alpha = 0.3) +
  geom_rect(data = geneExons, mapping = aes(xmin = start, xmax = end, ymin = 0.5, ymax = 0.7), position = "identity", stat = "identity", fill= "blue", alpha = 0.3) +
  geom_rect(data = geneExons, mapping = aes(xmin = start, xmax = end, ymin = -0.2, ymax = 0), position = "identity", stat = "identity", fill= "blue", alpha = 0.3) +
  geom_rect(data = geneProteinDomain, mapping = aes(xmin = start, xmax = end, ymin = 0, ymax = 0.5), position = "identity", stat = "identity", fill= "red", alpha = 0.8) +
  geom_text(data = geneProteinDomain, mapping = aes(x = start, y = 0.5, label = name, hjust = 0, vjust = 1)) +
  scale_x_continuous(name = "Exon Rank", breaks = geneExons$start, labels = geneExons$rank) +
  theme(axis.ticks = element_blank(), panel.background = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank()) +
  ggtitle(geneName)
  

####### PROTEIN DOMAIN COLOURS
proteinDomain = read.table(file = "~/hmf/analysis/fusions/SVA_VIS_PROTEIN_DOMAINS.tsv", sep = "\t", header = T) %>% 
  select(Transcript, Info) %>% distinct() %>%
  group_by(Info) %>% count() %>% arrange(-n) %>% select(Info)
proteinDomain$Info[1:10]


####### FUSIONS

fusions = read.table("/Users/jon/hmf/analysis/CPCT02390005T/CPCT02390005T.linx.fusions.tsv", sep = '\t', header  =T)

fusionDetails = read.table("/Users/jon/hmf/analysis/CPCT02390005T/CPCT02390005T.linx.fusions_detailed.csv", sep = ',', header  =T)

fusionDetailsSimple = fusionDetails %>%
  select(SampleId, ClusterId, GeneNameUp, ChrUp, PosUp, RegionTypeUp, StrandUp, ExonUp, ChrDown, PosDown, RegionTypeDown, StrandDown, ExonDown)



####### PLOTTING LEGEND


proteinDomain = read.table("/Users/jon/hmf/analysis/fusions/data/CPCT02050396T.cluster012.COMPLEX.sv374.protein_domain.circos") %>% 
  separate(V5, c("r", "g","b","alpha", "name"),sep=",")%>%
  mutate(
    chromosome = substring(V1, 3),
    r = substring(r, 13),
    alpha = substring(r, 1, nchar(alpha) - 1),
    name = substring(name, 6),
    name = gsub("[.]"," ", name),
    color = rgb(r,g,b,maxColorValue = 255)
  ) %>% 
  arrange(nchar(name)) 

proteinDomain$name = factor(proteinDomain$name, levels = unique(proteinDomain$name))


proteinDomainColors = proteinDomain %>% distinct(name, color)
proteinDomainColorsValues = setNames(proteinDomainColors$color, proteinDomainColors$name)

p1 = ggplot(proteinDomain) +
  geom_bar(aes(V1, V4, fill = name), stat = "identity") +
  scale_fill_manual(name = "Protein Domains", values = proteinDomainColorsValues) +
  guides(fill = guide_legend(ncol = 1)) +
  theme(legend.position = "left", legend.title = element_text(size=6), legend.text = element_text(size=5), 
        legend.key.size = unit(4,"mm"), 
        #legend.key.width = unit(4,"mm"),
        legend.spacing.x = unit(0,"mm"), legend.spacing.y = unit(0,"mm"))

legend <- cowplot::get_legend(p1)

ggsave("/Users/jon/hmf/analysis/fusions/testLegend.png", legend, height = 25*4, units = "mm")


ggdraw() + draw_plot(legend)


p <- ggplot(iris, aes(x=Sepal.Length, fill=Species)) + geom_density(alpha = 0.7)
ggdraw() +
  draw_image("/Users/jon/hmf/analysis/fusions/plot/CPCT02050396T.cluster012.COMPLEX.sv374.png") +
  draw_image("/Users/jon/hmf/analysis/fusions/testLegend.png")

