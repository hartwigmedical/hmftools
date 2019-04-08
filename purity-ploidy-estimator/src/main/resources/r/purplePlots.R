library(ggplot2)
library(dplyr)
theme_set(theme_bw())


# Parse the arguments
args <- commandArgs(trailing=T)
sample <- args[1]
purpleDir <- args[2]
plotDir   <- args[3]


purityPloidyRangePlot <- function(range) {

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
        scale_fill_gradientn(colours=rev(rainbow(1000, start=0, end=0.22))) +
        geom_segment(aes(y = 0.085, yend = 1.05, x=bestPloidy, xend = bestPloidy), linetype = "dashed") +
        geom_label(data = data.frame(), aes(x = bestPloidy, y = 1.05, label = round(bestPloidy, 2)), size = 2.5) +

        geom_segment(aes(y = bestPurity, yend = bestPurity, x=minPloidy, xend = maxPloidy + 0.4), linetype = "dashed") +
        geom_label(data = data.frame(), aes(y = bestPurity, x = maxPloidy + 0.4, label = paste0(bestPurity*100,"%" )), size = 2.5, hjust = 0.7) +
        theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title=element_text(size=6), legend.text=element_text(size=6)) +
        scale_y_continuous(labels = c("25%", "50%", "75%", "100%"), breaks = c(0.25, 0.5, 0.75, 1)) +
        xlab("Ploidy") + ylab("Purity") +  ggtitle("Purity/Ploidy Scores")

    return (result)
}

fittedSegmentsPlot <- function(fittedSegments) {
    fittedSegments = fittedSegments %>%
        filter(germlineStatus == "DIPLOID", bafCount > 0) %>%
        arrange(majorAllelePloidy) %>%
        mutate(
        Score = deviationPenalty * eventPenalty,
        Weight = bafCount,
        WeightedMajorAllelePloidyCumSum = cumsum(Weight * majorAllelePloidy),
        WeightedMajorAllelePloidyCumSumProportion = WeightedMajorAllelePloidyCumSum / max(WeightedMajorAllelePloidyCumSum))

    maxMajorAllelePloidy = ceiling(max(fittedSegments %>% filter(WeightedMajorAllelePloidyCumSumProportion <= 0.9) %>% select(majorAllelePloidy)))
    maxMinorAllelePloidy = maxMajorAllelePloidy - 1

    p = ggplot(fittedSegments, aes(x=majorAllelePloidy,y=minorAllelePloidy)) +
        geom_point(aes(size = Weight, color = Score), alpha = 0.7) +
        xlab("Major Allele") + ylab("Minor Allele") + ggtitle("Fitted Segment Scores") +
        scale_x_continuous(breaks = c(-200:200), limits = c(-0.1, maxMajorAllelePloidy)) +
        scale_y_continuous(breaks = c(-200:200), limits = c(-0.1, maxMinorAllelePloidy)) +
        scale_color_gradientn(colours=rev(rainbow(1000, start=0, end=0.25))) +
        theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), legend.position = "right", legend.title=element_text(size=6), legend.text=element_text(size=6)) +
        scale_size(range = c(1,9), guide = "none")

    return (p)

}

rangeDF = read.table(file = paste0(purpleDir, "/", sample, ".purple.purity.range"), sep = "\t", header = T, comment.char = "!") %>% select(Purity = X.Purity, Ploidy, Score)
rangePlot = purityPloidyRangePlot(rangeDF)
ggsave(filename = paste0(plotDir, "/", sample, ".purity.range.png"), rangePlot, units = "in", height = 4, width = 4.8, scale = 1)


fittedSegmentsDF = read.table(file = paste0(purpleDir, "/", sample, ".purple.fitted"), sep = "\t", header = T, comment.char = "!")
fittedSegmentsPlot = fittedSegmentsPlot(fittedSegmentsDF)
ggsave(filename = paste0(plotDir, "/", sample, ".fitted.segments.png"), fittedSegmentsPlot, units = "in", height = 4, width = 4.8, scale = 1)
