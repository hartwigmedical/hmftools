library(ggplot2)
library(dplyr)
theme_set(theme_bw())

# Parse the arguments
args <- commandArgs(trailing=T)
sample <- args[1]
purpleDir <- args[2]
plotDir   <- args[3]

purity_ploidy_range_plot <- function(range) {

    bestPurity = range[1, "Purity"]
    bestPloidy = range[1, "Ploidy"]

    range =  range %>%
        arrange(Purity, Ploidy) %>%
        group_by(Purity) %>%
        mutate(
        Score = pmin(4, Score),
        LeftPloidy = lag(Ploidy),
        RightPloidy = lead(Ploidy),
        xmin = Ploidy - (Ploidy - LeftPloidy) / 2,
        xmax = Ploidy + (RightPloidy - Ploidy) / 2,
        ymin = Purity - 0.005,
        ymax = Purity + 0.005,
        xmin = ifelse(is.na(xmin), Ploidy, xmin),
        xmax = ifelse(is.na(xmax), Ploidy, xmax))

    maxPloidy = min(range %>% arrange(Purity, -Ploidy) %>% group_by(Purity)  %>% filter(row_number() == 1) %>% select(Purity, Ploidy = xmax) %>% ungroup() %>% select(Ploidy))
    minPloidy = max(range %>% arrange(Purity, Ploidy) %>% group_by(Purity)  %>% filter(row_number() == 1) %>% select(Purity, MaxPloidy = xmin) %>% ungroup() %>% select(MaxPloidy))

    range = range %>%
        filter(xmin <= maxPloidy, xmax >= minPloidy) %>%
        mutate(xmax = pmin(xmax, maxPloidy), xmin = pmax(xmin, minPloidy))

    result = ggplot(range) +
        geom_rect(aes(fill=Score, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
        scale_fill_gradientn(colours=c("blue","green","yellow","orange", "red"), limits = c(min(range$Score), 4)) +
        geom_segment(aes(y = 0.085, yend = 1.05, x=bestPloidy, xend = bestPloidy), linetype = "dashed") +
        geom_label(data = data.frame(), aes(x = bestPloidy, y = 1.05, label = round(bestPloidy, 2)), size = 2.5) +

        geom_segment(aes(y = bestPurity, yend = bestPurity, x=minPloidy, xend = maxPloidy + 0.4), linetype = "dashed") +
        geom_label(data = data.frame(), aes(y = bestPurity, x = maxPloidy + 0.4, label = paste0(bestPurity*100,"%" )), size = 2.5, hjust = 0.7) +
        theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title=element_text(size=6), legend.text=element_text(size=6)) +
        scale_y_continuous(labels = c("25%", "50%", "75%", "100%"), breaks = c(0.25, 0.5, 0.75, 1)) +
        xlab("Ploidy") + ylab("Purity") +  ggtitle("Purity/Ploidy Scores")

    return (result)
}

fitted_segments_plot <- function(fittedSegments) {
    fittedSegments = fittedSegments %>%
        filter(germlineStatus == "DIPLOID", bafCount > 0) %>%
        arrange(majorAllelePloidy) %>%
        mutate(
        Score = deviationPenalty * eventPenalty,
        Weight = bafCount,
        WeightedMajorAllelePloidyCumSum = cumsum(Weight * majorAllelePloidy),
        WeightedMajorAllelePloidyCumSumProportion = WeightedMajorAllelePloidyCumSum / max(WeightedMajorAllelePloidyCumSum))

    maxData = fittedSegments %>% filter(WeightedMajorAllelePloidyCumSumProportion <= 0.9) %>% select(majorAllelePloidy, Score)
    maxScore = ceiling(max(maxData$Score))
    minScore = floor(min(maxData$Score))
    maxMajorAllelePloidy = ceiling(max(maxData$majorAllelePloidy))
    maxMinorAllelePloidy = maxMajorAllelePloidy - 1

    p = ggplot(fittedSegments, aes(x=majorAllelePloidy,y=minorAllelePloidy)) +
        geom_point(aes(size = Weight, color = Score), alpha = 0.7) +
        xlab("Major Allele") + ylab("Minor Allele") + ggtitle("Fitted Segment Scores") +
        scale_x_continuous(breaks = c(-200:200), limits = c(-0.1, maxMajorAllelePloidy)) +
        scale_y_continuous(breaks = c(-200:200), limits = c(-0.1, maxMinorAllelePloidy)) +
        scale_color_gradientn(colours=c("blue","green","yellow","orange", "red"), limits = c(minScore, maxScore)) +
        theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title=element_text(size=6), legend.text=element_text(size=6)) +
        scale_size(range = c(1,9), guide = "none")

    return (p)

}



minor_allele_ploidy_pdf <- function(copyNumberRegions) {
    cnColours = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69")
    cnColours = setNames(cnColours, c("CN0", "CN1","CN2","CN3","CN4", "CN5", "CN6+"))

    totalBafCount = sum(copyNumberRegions$bafCount)

    maxAllelePloidy = copyNumberRegions %>%
        mutate(bucket = ceiling(minorAllelePloidy)) %>%
        group_by(bucket) %>%
        summarise(n = sum(bafCount)) %>%
        mutate(cumn = cumsum(n), proportion =  cumn / max(cumn)) %>%
        arrange(proportion) %>%
        filter(proportion > 0.9) %>%
        filter(row_number() == 1) %>%
        pull(bucket)

    copyNumberRegions = copyNumberRegions %>%
        filter(!chromosome %in% c('X','Y'), bafCount > 0) %>%
        mutate(
        cn = round(copyNumber),
        cn = ifelse(cn >=6, "CN6+", paste0("CN", cn)),
        weight = bafCount/totalBafCount )

    ggplot(copyNumberRegions, aes(x = minorAllelePloidy)) +
        geom_histogram(aes(weight = bafCount, fill = cn), alpha =1, binwidth = 0.1, color = "black",  position = "stack") +
        scale_x_continuous(breaks = c(0:10), limits = c(-0.1, maxAllelePloidy + 0.1)) +
        scale_fill_manual(values = cnColours) +
        theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title = element_blank()) +
        xlab("Minor Allele Ploidy") + ylab("Baf Count") + ggtitle("Minor Allele Ploidy PDF")
}


copynumber_pdf <- function(copyNumberRegions) {

    mapColours = c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462")
    mapColours = setNames(mapColours, c("MAP0", "MAP1","MAP2","MAP3","MAP4", "MAP5+"))

    totalBafCount = sum(copyNumberRegions$bafCount)

    maxCopyNumber = copyNumberRegions %>%
        mutate(bucket = ceiling(copyNumber)) %>%
        group_by(bucket) %>%
        summarise(n = sum(bafCount)) %>%
        mutate(cumn = cumsum(n), proportion =  cumn / max(cumn)) %>%
        arrange(proportion) %>%
        filter(proportion > 0.9) %>%
        filter(row_number() == 1) %>%
        pull(bucket)

    copyNumberRegions = copyNumberRegions %>%
        filter(!chromosome %in% c('X','Y'), bafCount > 0) %>%
        mutate(
        map = round(minorAllelePloidy),
        map = ifelse(map>=5, "MAP5+", paste0("MAP", map)),
        map = paste0("MAP", round(minorAllelePloidy)),
        chromosome = factor(chromosome, levels= c(1:22), ordered = T),
        weight = bafCount/totalBafCount )

    ggplot(copyNumberRegions, aes(x = copyNumber)) +
        geom_histogram(aes(weight = bafCount, fill = map), alpha = 1,  binwidth = 0.1, color = "black",  position = "stack") +
        scale_fill_manual(values = mapColours) +
        scale_x_continuous(breaks = c(0:10), limits = c(-0.1, maxCopyNumber + 0.1)) +
        theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.title = element_blank()) +
        xlab("Copy Number") + ylab("Baf Count") + ggtitle("Copy Number PDF")
}

copyNumbers = read.table(file = paste0(purpleDir, "/", sample, ".purple.cnv"), sep = "\t", header = T, comment.char = "!") %>% select(chromosome = X.chromosome, everything())
copyNumberPDF = copynumber_pdf(copyNumbers)
ggsave(filename = paste0(plotDir, "/", sample, ".copynumber.png"), copyNumberPDF, units = "in", height = 4, width = 4.8, scale = 1)

minorAllelePloidyPDF = minor_allele_ploidy_pdf(copyNumbers)
ggsave(filename = paste0(plotDir, "/", sample, ".map.png"), minorAllelePloidyPDF, units = "in", height = 4, width = 4.8, scale = 1)

rangeDF = read.table(file = paste0(purpleDir, "/", sample, ".purple.purity.range"), sep = "\t", header = T, comment.char = "!") %>% select(Purity = X.Purity, Ploidy, Score)
rangePlot = purity_ploidy_range_plot(rangeDF)
ggsave(filename = paste0(plotDir, "/", sample, ".purity.range.png"), rangePlot, units = "in", height = 4, width = 4.8, scale = 1)


fittedSegmentsDF = read.table(file = paste0(purpleDir, "/", sample, ".purple.fitted"), sep = "\t", header = T, comment.char = "!")
fittedSegmentsPlot = fitted_segments_plot(fittedSegmentsDF)
ggsave(filename = paste0(plotDir, "/", sample, ".fitted.segments.png"), fittedSegmentsPlot, units = "in", height = 4, width = 4.8, scale = 1)
