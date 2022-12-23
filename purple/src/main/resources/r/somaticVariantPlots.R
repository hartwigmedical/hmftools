library(VariantAnnotation)
library(ggplot2)
library(dplyr)
theme_set(theme_bw())

#sample = "COLO829v003T"
#purpleDir <- "~/hmf/analysis/COLO829T/purple"
#plotDir   <- "~/hmf/analysis/COLO829T/purple/plot"

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
    ref = as.character(ref(vcf)), 
    alt = as.character(vcf.alt),  
    filter = as.character(vcf.rowRanges$FILTER), 
    minorAllelePloidy = vcf.info$PURPLE_MACN,
    ploidy = vcf.info$PURPLE_VCN,
    copyNumber = vcf.info$PURPLE_CN,
    kataegis = vcf.info$KT, stringsAsFactors = F)
  
  return (vcf.df)
}

standard_mutation <- function(types) {
  types = gsub("G>T", "C>A", types)
  types = gsub("G>C", "C>G", types)
  types = gsub("G>A", "C>T", types)
  types = gsub("A>T", "T>A", types)
  types = gsub("A>G", "T>C", types)
  types = gsub("A>C", "T>G", types)
  return(types)
}

clonality_plot <- function(somaticBuckets, clonalityModel) {
  clonalityVariants = somaticBuckets %>% group_by(variantCopyNumberBucket) %>% summarise(count = sum(count))
  
  subclonalPercentage = clonalityModel %>% 
    group_by(bucket) %>% 
    mutate(totalWeight = sum(bucketWeight)) %>% 
    filter(isSubclonal) %>%
    summarise(
      isSubclonal = T,
      bucketWeight = sum(bucketWeight), 
      subclonalLikelihood = ifelse(bucketWeight == 0, 0, bucketWeight / max(totalWeight)))

  nonResidualModel = clonalityModel %>% filter(peak != 0)

  nonResidualSubclonalPercentage = nonResidualModel %>%
    group_by(bucket) %>%
    mutate(totalWeight = sum(bucketWeight)) %>%
    filter(isSubclonal) %>%
    summarise(
    isSubclonal = T,
    bucketWeight = sum(bucketWeight),
    subclonalLikelihood = ifelse(bucketWeight == 0, 0, bucketWeight / max(totalWeight)))

  combinedModel = nonResidualModel %>%
    group_by(bucket) %>% 
    summarise(bucketWeight = sum(bucketWeight))
  
  singleBlue = "#6baed6"
  singleRed = "#d94701"
  
  pTop = ggplot() +
    geom_bar(data=clonalityVariants, aes(x = variantCopyNumberBucket, weight = count), fill=singleBlue, col=singleBlue,  alpha = .4, size = 0.07, width = 0.05) +
    geom_line(data=combinedModel , aes(x = bucket, y = bucketWeight), position = "identity", alpha = 0.8) +
    geom_line(data=nonResidualModel, aes(x = bucket, y = bucketWeight, color = peak), position = "identity") +
    geom_area(data=nonResidualSubclonalPercentage %>% filter(isSubclonal), aes(x = bucket, y = bucketWeight), position = "identity",  alpha = 0.3, fill = singleRed, color = singleRed) +
    ggtitle("") + xlab("Variant Copy Number") + ylab("") +
    scale_y_continuous(expand=c(0.02, 0.02)) +
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position="none") +
    scale_x_continuous( expand=c(0.01, 0.01), limits = c(0, 3.5)) 

  pBottom = ggplot(data = subclonalPercentage) +
    geom_bar(width = 0.05, aes(x = bucket, y = subclonalLikelihood), stat = "identity", fill=singleRed, col=singleRed,  alpha = 0.3) + 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank(), axis.ticks = element_blank()) +
    xlab("") + ylab("") +
    scale_y_continuous(labels = c("0%", "25%","50%","75%","100%"), breaks = c(0, 0.25, 0.5, 0.75, 1), expand=c(0.02, 0.02), limits = c(0, 1)) +
    scale_x_continuous( expand=c(0.01, 0.01), limits = c(0, 3.5)) 
  
  pFinal = cowplot::plot_grid(pTop, pBottom, ncol = 1, rel_heights = c(5, 1), align = "v")
  return(pFinal)
}

maxRainfallSNVs=100000

rainfall_plot <- function(somaticVariants) {
  
  strandColours = c("#6bd692", "#7e6bd6")
  strandColours = setNames(strandColours, c("Forward", "Reverse"))
  
  singleSubstitutionColours = c("#14B0EF","#060809","#E00714","#BFBEBF","#90CA4B","#E9BBB8")
  singleSubstitutionColours = setNames(singleSubstitutionColours, c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
  
  snps = somaticVariants %>% 
    filter(filter == 'PASS', nchar(ref) == 1, nchar(alt) == 1) %>% 
    mutate(
      mutation = paste(ref, alt, sep = ">"), 
      mutation = standard_mutation(mutation),
      prevPos = lag(pos), 
      nextPos = lead(pos),
      prevPos = ifelse(is.na(prevPos), 0, prevPos),
      nextPos = ifelse(is.na(nextPos), 0, nextPos),
      prevDis = abs(pos - prevPos),
      nextDis = abs(nextPos - pos),
      distanceToNeighbour = pmin(prevDis, nextDis),
      rank = row_number())

  if(nrow(snps) > maxRainfallSNVs) {
    snps = head(snps,maxRainfallSNVs)
  }
  
  kataegis = snps %>% mutate(ymin = min(distanceToNeighbour), ymax = max(distanceToNeighbour)) %>% 
    filter(!is.na(kataegis)) %>% 
    group_by(kataegis, ymin, ymax) %>%
    summarise(xmin = min(rank), xmax = max(rank)) %>%
    mutate(strand = ifelse(substring(kataegis, 1, 3) == "FWD", "Forward", "Reverse"))

  p = ggplot() +
    geom_rect(data = kataegis, mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = strand), alpha = 0.6) + 
    geom_point(data = snps, mapping = aes(x = rank, y = distanceToNeighbour, color = mutation), size = 0.1) +
    scale_y_log10(labels = function(x) {format(x, scientific = FALSE)}) + 
    ylab("Intermutation distance (bp)") + xlab("Mutation number") +
    scale_color_manual(values = singleSubstitutionColours, name = "Mutation") + 
    scale_fill_manual(values = strandColours, name = "Kataegis Regions") + 
    guides(color = guide_legend(override.aes = list(size = 2)))
  
  return (p)
}

somatic_ploidy_pdf <- function(somaticBuckets) {
  cnColours = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69")
  cnColours = setNames(cnColours, c("CN0", "CN1","CN2","CN3","CN4", "CN5", "CN6+"))
  
  maxPloidy = somaticBuckets %>%
    mutate(bucket = ceiling(variantCopyNumberBucket)) %>%
    group_by(bucket) %>%
    summarise(n = sum(count)) %>%
    mutate(cumn = cumsum(n), proportion =  cumn / max(cumn)) %>%
    arrange(proportion) %>%
    filter(proportion > 0.95)   %>%
    filter(row_number() == 1) %>%
    pull(bucket)
  
  somatics = somaticBuckets %>%
    mutate(
      cn = pmax(0, round(copyNumberBucket)),
      cn = ifelse(cn >=6, "CN6+", paste0("CN", cn)))
  
  ggplot(somatics, aes(x = variantCopyNumberBucket)) +
    geom_bar(aes(fill = cn, weight = count), alpha=1, color = "black",  position = "stack", size = 0.07, width = 0.05) +
    scale_x_continuous(breaks = c(0:10), limits = c(-0.1, maxPloidy + 1.1)) +
    scale_fill_manual(values = cnColours) +
    theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title = element_blank()) +
    xlab("Variant Copy Number") + ylab("Count") + ggtitle("Somatic Variant Copy Number PDF")
}

somaticBuckets = read.table(paste0(purpleDir, "/", sample, ".purple.somatic.hist.tsv"), sep = "\t", header = T, numerals = "no.loss", skipNul = T)
somaticVariantPDF = somatic_ploidy_pdf(somaticBuckets)
ggsave(filename = paste0(plotDir, "/", sample, ".somatic.png"), somaticVariantPDF, units = "in", height = 4, width = 4.8, scale = 1)

clonalityModel = read.table(paste0(purpleDir, "/", sample, ".purple.somatic.clonality.tsv"), sep = "\t", header = T, numerals = "no.loss", skipNul = T) %>%
  mutate(isSubclonal = isSubclonal == "true", isValid = isValid == "true", peak = as.character(peak), bucketWeight = as.numeric(as.character(bucketWeight))) %>% filter(isValid)
clonalityModelPlot = clonality_plot(somaticBuckets, clonalityModel)
ggsave(filename = paste0(plotDir, "/", sample, ".somatic.clonality.png"), clonalityModelPlot, units = "in", height = 6, width = 8, scale = 1)

vcf = readVcf(paste0(purpleDir, "/", sample, ".purple.somatic.vcf.gz"))
somaticVariants = vcf_data_frame(vcf)
rainfallPlot = rainfall_plot(somaticVariants)
ggsave(filename = paste0(plotDir, "/", sample, ".somatic.rainfall.png"), rainfallPlot, units = "in", height = 4, width = 8, scale = 1)
