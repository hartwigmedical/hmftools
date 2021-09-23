library(ggplot2)
library(dplyr)
theme_set(theme_bw())

#sample = "COLO829v003T"
#purpleDir = "~/hmf/analysis/COLO829T/purple"
#plotDir = "~/hmf/analysis/COLO829T/purple/plot"

# Parse the arguments
args <- commandArgs(trailing=T)
sample <- args[1]
purpleDir <- args[2]
plotDir   <- args[3]

purity_ploidy_range_plot <- function(bestFit, range) {

    bestPurity = bestFit[1, "purity"]
    bestPloidy = bestFit[1, "ploidy"]
    bestScore = bestFit[1, "score"]
    
    range =  range %>%
        arrange(purity, ploidy) %>%
        group_by(purity) %>%
        mutate(
        absScore = pmin(4, score),
        score = pmin(1, abs(score - bestScore) / score),
        leftPloidy = lag(ploidy),
        rightPloidy = lead(ploidy),
        xmin = ploidy - (ploidy - leftPloidy) / 2,
        xmax = ploidy + (rightPloidy - ploidy) / 2,
        ymin = purity - 0.005,
        ymax = purity + 0.005,
        xmin = ifelse(is.na(xmin), ploidy, xmin),
        xmax = ifelse(is.na(xmax), ploidy, xmax))

    maxPloidy = min(range %>% arrange(purity, -ploidy) %>% group_by(purity)  %>% filter(row_number() == 1) %>% select(purity, ploidy = xmax) %>% ungroup() %>% select(ploidy))
    minPloidy = max(range %>% arrange(purity, ploidy) %>% group_by(purity)  %>% filter(row_number() == 1) %>% select(purity, maxPloidy = xmin) %>% ungroup() %>% select(maxPloidy))

    maxPloidy = max(maxPloidy, bestPloidy)
    minPloidy = min(minPloidy, bestPloidy)
    
    range = range %>%
        filter(xmin <= maxPloidy, xmax >= minPloidy) %>%
        mutate(xmax = pmin(xmax, maxPloidy), xmin = pmax(xmin, minPloidy))

    result = ggplot(range) +
        geom_rect(aes(fill=score, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
        scale_fill_gradientn(colours=c("blue","blue", "green", "yellow","orange", "red", "red2"), limits = c(0, 1), values=c(0, 0.0999, 0.1, 0.5, 0.8, 0.9, 1), breaks = c(0.1, 0.5, 1), labels = c("10%", "50%", "100%"), name = "Relative\nScore") +
        geom_segment(aes(y = 0.085, yend = 1.05, x=bestPloidy, xend = bestPloidy), linetype = "dashed", size = 0.1) +
        geom_label(data = data.frame(), aes(x = bestPloidy, y = 1.05, label = round(bestPloidy, 2)), size = 2.5) +
        geom_segment(aes(y = bestPurity, yend = bestPurity, x=minPloidy, xend = maxPloidy + 0.4), linetype = "dashed", size = 0.1) +
        geom_label(data = data.frame(), aes(y = bestPurity, x = maxPloidy + 0.4, label = paste0(bestPurity*100,"%" )), size = 2.5, hjust = 0.7) +
        theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title=element_text(size=6), legend.text=element_text(size=6)) +
        scale_y_continuous(labels = c("25%", "50%", "75%", "100%"), breaks = c(0.25, 0.5, 0.75, 1)) +
        xlab("Ploidy") + ylab("Purity") +  ggtitle("Purity/Ploidy Scores")

    return (result)
}

fitted_segments_plot <- function(fittedSegments) {
    fittedSegments = fittedSegments %>%
        filter(germlineStatus == "DIPLOID", bafCount > 0) %>%
        arrange(majorAlleleCopyNumber) %>%
        mutate(
        Score = deviationPenalty * eventPenalty,
        Weight = bafCount,
        WeightedMajorAllelePloidyCumSum = cumsum(Weight * majorAlleleCopyNumber),
        WeightedMajorAllelePloidyCumSumProportion = WeightedMajorAllelePloidyCumSum / max(WeightedMajorAllelePloidyCumSum))

    maxData = fittedSegments %>% filter(WeightedMajorAllelePloidyCumSumProportion <= 0.9) %>% select(majorAlleleCopyNumber, Score)
    maxScore = ceiling(max(maxData$Score))
    minScore = floor(min(maxData$Score))
    minMajorAllelePloidy = min(0, floor(min(maxData$majorAlleleCopyNumber)))
    maxMajorAllelePloidy = ceiling(max(maxData$majorAlleleCopyNumber))
    maxMinorAllelePloidy = maxMajorAllelePloidy - 1
    
    p = ggplot(fittedSegments, aes(x=majorAlleleCopyNumber,y=minorAlleleCopyNumber)) +
        geom_point(aes(size = Weight, color = Score), alpha = 0.7) +
        xlab("Major Allele") + ylab("Minor Allele") + ggtitle("Segment Scores") +
        scale_x_continuous(breaks = c(-200:200), limits = c(minMajorAllelePloidy, maxMajorAllelePloidy)) +
        scale_y_continuous(breaks = c(-200:200), limits = c(0, maxMinorAllelePloidy)) +
        scale_color_gradientn(colours=c("blue","green","yellow","orange", "red"), limits = c(minScore, maxScore)) +
        theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title=element_text(size=6), legend.text=element_text(size=6)) +
        scale_size(range = c(1,9), guide = "none")

    return (p)

}

minor_allele_ploidy_pdf <- function(copyNumberRegions) {
    cnColours = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69")
    cnColours = setNames(cnColours, c("CN0", "CN1","CN2","CN3","CN4", "CN5", "CN6+"))

    totalBafCount = sum(copyNumberRegions$bafCount)

    copyNumberRegions = copyNumberRegions %>%
      filter(!chromosome %in% c('X','Y'), bafCount > 0) %>%
      mutate(
        cn = pmax(0, round(copyNumber)),
        cn = ifelse(cn >=6, "CN6+", paste0("CN", cn)),
        weight = bafCount/totalBafCount )
    
    maxAllelePloidy = copyNumberRegions %>%
        mutate(bucket = ceiling(minorAlleleCopyNumber)) %>%
        group_by(bucket) %>%
        summarise(n = sum(bafCount)) %>%
        mutate(cumn = cumsum(n), proportion =  cumn / max(cumn)) %>%
        arrange(proportion) %>%
        filter(proportion > 0.9) %>%
        filter(row_number() == 1) %>%
        pull(bucket)

    ggplot(copyNumberRegions, aes(x = minorAlleleCopyNumber)) +
        geom_histogram(aes(weight = bafCount, fill = cn), alpha =1, binwidth = 0.1, color = "black",  position = "stack", size = 0.07) +
        scale_x_continuous(breaks = c(0:10), limits = c(-0.1, maxAllelePloidy + 0.1)) +
        scale_fill_manual(values = cnColours) +
        theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title = element_blank()) +
        xlab("Minor Allele Copy Number") + ylab("Baf Count") + ggtitle("Minor Allele Copy Number PDF")
}

copynumber_pdf <- function(copyNumberRegions) {

    mapColours = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462")
    mapColours = setNames(mapColours, c("MACN0", "MACN1","MACN2","MACN3","MACN4", "MACN5+"))

    totalBafCount = sum(copyNumberRegions$bafCount)

    copyNumberRegions = copyNumberRegions %>%
      filter(!chromosome %in% c('X','Y'), bafCount > 0) %>%
      mutate(
        map = round(minorAlleleCopyNumber),
        map = ifelse(map>=5, "MACN5+", paste0("MACN", map)),
        chromosome = factor(chromosome, levels= c(1:22), ordered = T),
        weight = bafCount/totalBafCount )
    
    maxCopyNumber = copyNumberRegions %>%
        mutate(bucket = ceiling(copyNumber)) %>%
        group_by(bucket) %>%
        summarise(n = sum(bafCount)) %>%
        mutate(cumn = cumsum(n), proportion =  cumn / max(cumn)) %>%
        arrange(proportion) %>%
        filter(proportion > 0.9) %>%
        filter(row_number() == 1) %>%
        pull(bucket)
    
    minCopyNumber = floor(min(copyNumberRegions$copyNumber))
    
    ggplot(copyNumberRegions, aes(x = copyNumber)) +
        geom_histogram(aes(weight = bafCount, fill = map), alpha = 1,  binwidth = 0.1, color = "black",  position = "stack", size = 0.07) +
        scale_fill_manual(values = mapColours) +
        scale_x_continuous(breaks = c((minCopyNumber - 1):(maxCopyNumber + 1)), limits = c(minCopyNumber - 0.1, maxCopyNumber + 0.1)) +
        theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.title = element_blank()) +
        xlab("Copy Number") + ylab("Baf Count") + ggtitle("Copy Number PDF")
}

copyNumbers = read.table(file = paste0(purpleDir, "/", sample, ".purple.cnv.somatic.tsv"), sep = "\t", header = T, comment.char = "!") %>%
  mutate(chromosome = gsub("chr", "", chromosome)) %>%
  filter(!chromosome %in% c('X','Y'), bafCount > 0) 

if (nrow(copyNumbers) > 0) {
  copyNumberPDF = copynumber_pdf(copyNumbers)
  ggsave(filename = paste0(plotDir, "/", sample, ".copynumber.png"), copyNumberPDF, units = "in", height = 4, width = 4.8, scale = 1)
  
  minorAllelePloidyPDF = minor_allele_ploidy_pdf(copyNumbers)
  ggsave(filename = paste0(plotDir, "/", sample, ".map.png"), minorAllelePloidyPDF, units = "in", height = 4, width = 4.8, scale = 1)  
}


bestFitDF = read.table(file = paste0(purpleDir, "/", sample, ".purple.purity.tsv"), sep = "\t", header = T, comment.char = "!") %>% select(purity, ploidy, score)
rangeDF = read.table(file = paste0(purpleDir, "/", sample, ".purple.purity.range.tsv"), sep = "\t", header = T, comment.char = "!") %>%
    select(purity, ploidy, score)
rangePlot = purity_ploidy_range_plot(bestFitDF, rangeDF)
ggsave(filename = paste0(plotDir, "/", sample, ".purity.range.png"), rangePlot, units = "in", height = 4, width = 4.8, scale = 1)

fittedSegmentsDF = read.table(file = paste0(purpleDir, "/", sample, ".purple.segment.tsv"), sep = "\t", header = T, comment.char = "!")
fittedSegmentsPlot = fitted_segments_plot(fittedSegmentsDF)
ggsave(filename = paste0(plotDir, "/", sample, ".segment.png"), fittedSegmentsPlot, units = "in", height = 4, width = 4.8, scale = 1)

