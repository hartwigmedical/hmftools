library(ggplot2)
library(dplyr)
library(VariantAnnotation)
theme_set(theme_bw())


# Parse the arguments
args <- commandArgs(trailing=T)
sample <- args[1]
purpleDir <- args[2]
plotDir   <- args[3]


vcf_data_frame<- function(vcf) {
  vcf.rowRanges = rowRanges(vcf)
  vcf.info = info(vcf)
  vcf.alt = CharacterList(alt(vcf))
  vcf.alt[elementNROWS(vcf.alt) > 1L ] <- lapply(vcf.alt[elementNROWS(vcf.alt) > 1L ], paste0, collapse=",")
  
  vcf.df = data.frame(
    chromosome = seqnames(vcf), 
    pos = start(vcf), 
    ref = ref(vcf), 
    alt = as.character(vcf.alt),  
    filter = as.character(vcf.rowRanges$FILTER), 
    minorAllelePloidy = vcf.info$PURPLE_MAP, 
    ploidy = vcf.info$PURPLE_PLOIDY, 
    copyNumber = vcf.info$PURPLE_CN)
  
  return (vcf.df)
}

somatic_ploidy_pdf <- function(somatics) {
  cnColours = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69")
  cnColours = setNames(cnColours, c("CN0", "CN1","CN2","CN3","CN4", "CN5", "CN6+"))
  
  maxPloidy = somatics %>%
    mutate(bucket = ceiling(ploidy)) %>%
    group_by(bucket) %>%
    summarise(n = n()) %>%
    mutate(cumn = cumsum(n), proportion =  cumn / max(cumn)) %>%
    arrange(proportion) %>%
    filter(proportion > 0.95) %>%
    filter(row_number() == 1) %>%
    pull(bucket)
  
  somatics = somatics %>%
    filter(!chromosome %in% c('MT'), filter == 'PASS') %>%
    mutate(
      cn = pmax(0, round(copyNumber)),
      cn = ifelse(cn >=6, "CN6+", paste0("CN", cn)))
  
  ggplot(somatics, aes(x = ploidy)) +
    geom_histogram(aes(fill = cn), alpha =1, binwidth = 0.1, color = "black",  position = "stack") +
    scale_x_continuous(breaks = c(0:10), limits = c(-0.1, maxPloidy + 1.1)) +
    scale_fill_manual(values = cnColours) +
    theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title = element_blank()) +
    xlab("Ploidy") + ylab("Count") + ggtitle("Somatic Variant Ploidy PDF")
}


somaticVariants = vcf_data_frame(readVcf(file = paste0(purpleDir, "/", sample, ".purple.somatic.vcf.gz")))
somaticVariantPDF = somatic_ploidy_pdf(somaticVariants)
ggsave(filename = paste0(plotDir, "/", sample, ".variant.png"), somaticVariantPDF, units = "in", height = 4, width = 4.8, scale = 1)
