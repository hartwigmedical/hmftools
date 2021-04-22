library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
theme_set(theme_bw())

ploidyPenaltyFactor = 0.4
standardDeviation = 0.05
minStandardDeviationPerPloidyPoint = 1.5
majorAlleleSubOnePenaltyMultiplier = 1
majorAlleleSubOneAdditionalPenalty = 1.5
baselineDeviation = 0.1

alleleDeviation <- function(purity, normFactor, ploidy) {
    ploidyDistanceFromInteger = 0.5
    if (ploidy > -0.5) {
        ploidyDistanceFromInteger = abs(ploidy - round(ploidy))
    }

    standardDeviationsPerPloidy = max(minStandardDeviationPerPloidyPoint, purity * normFactor / 2/ standardDeviation)
    return (2 * pnorm(ploidyDistanceFromInteger * standardDeviationsPerPloidy) - 1 + max(-0.5-ploidy,0))
}

subMininimumPloidyPenalty <- function(minPloidy, ploidy) {
    penalty = - majorAlleleSubOneAdditionalPenalty * (ploidy - minPloidy)
    return (min(majorAlleleSubOneAdditionalPenalty, max(penalty, 0)))
}

minorAlleleDeviation <- function(purity, normFactor, ploidy) {
    deviation = alleleDeviation(purity, normFactor, ploidy) + subMininimumPloidyPenalty(0, ploidy)
    return (max(deviation, baselineDeviation))
}

majorAlleleDeviation <- function(purity, normFactor, ploidy) {
    majorAlleleMultiplier = 1
    if (ploidy > 0 && ploidy < 1) {
        majorAlleleMultiplier = pmax(1,majorAlleleSubOnePenaltyMultiplier*(1-ploidy))
    }

    deviation = majorAlleleMultiplier * alleleDeviation(purity, normFactor, ploidy) + subMininimumPloidyPenalty(1, ploidy)
    return (max(deviation, baselineDeviation))
}

eventPenalty <- function(majorAllele, minorAllele) {
    wholeGenomeDoublingDistance = 1 + (abs(majorAllele - 2)) + (abs(minorAllele - 2));
    singleEventDistance = (abs(majorAllele - 1)) + (abs(minorAllele - 1));

    return (1 + ploidyPenaltyFactor * min(singleEventDistance, wholeGenomeDoublingDistance))
}

ploidy = seq(-1.2, 3, by = 0.01)
highPurityMajorPenalty = data.frame(ploidy, penalty = sapply(ploidy, function (x) {majorAlleleDeviation(0.7, 1, x) + 0.005}))
highPurityMajorPenalty$Penalty <- "Major"
highPurityMajorPenalty$purity <- 0.7

highPurityMinorPenalty = data.frame(ploidy, penalty = sapply(ploidy, function (x) {minorAlleleDeviation(0.7, 1, x) - 0.005}))
highPurityMinorPenalty$Penalty <- "Minor"
highPurityMinorPenalty$purity <- 0.7

lowPurityMajorPenalty = data.frame(ploidy, penalty = sapply(ploidy, function (x) {majorAlleleDeviation(0.3, 1, x) + 0.005}))
lowPurityMajorPenalty$Penalty <- "Major"
lowPurityMajorPenalty$purity <- 0.3

lowPurityMinorPenalty = data.frame(ploidy, penalty = sapply(ploidy, function (x) {minorAlleleDeviation(0.3, 1, x) - 0.005}))
lowPurityMinorPenalty$Penalty <- "Minor"
lowPurityMinorPenalty$purity <- 0.3

df =  bind_rows(highPurityMajorPenalty, highPurityMinorPenalty,lowPurityMajorPenalty, lowPurityMinorPenalty) %>%
mutate(Purity = factor(purity, levels = c(0.7, 0.3), ordered = T), Penalty = factor(Penalty, levels = c("Minor", "Major"), ordered = T))

ploidyColours = setNames(c("#6baed6", "#fb6a4a"),c("Major", "Minor"))


pDeviationPenalty = ggplot(df) +
    geom_line(aes(x = ploidy, y = penalty, color = Penalty, linetype = Purity)) +
    scale_color_manual(values = ploidyColours) +
    xlab("Ploidy") + ylab("Penalty") + ggtitle("Ploidy Deviation") +
    theme( panel.border = element_blank()) +
    theme(axis.ticks = element_blank(), legend.position="bottom") +
    ylim(0,3.2)
pDeviationPenalty

save_plot("~/hmf/repos/hmftools/purple/src/main/resources/readme/FittedPurityDeviationPenalty.png", pDeviationPenalty, base_width = 8, base_height = 4)


purityMatrix <- function(purity, ploidy) {
    resultMatrix = matrix(nrow = length(ploidy), ncol = length(ploidy))
    for (i in c(1:length(ploidy))){
        cat (i, " ")
        for (j in c(1:i)){
            majorPloidy = ploidy[i]
            minorPloidy = ploidy[j]
            totalPenalty = (majorAlleleDeviation(purity, 1, majorPloidy) + minorAlleleDeviation(purity, 1, minorPloidy)) * eventPenalty(majorPloidy, minorPloidy)
            resultMatrix[i,j] = totalPenalty
        }
    }

    return (resultMatrix)
}

purityDataFrame <- function(mat, ploidy) {
    df = cbind(majorAllele = ploidy, data.frame(mat))
    colnames(df) <- c("MajorAllele", ploidy)
    df = df %>% gather(MinorAllele, Penalty, -MajorAllele) %>%
        filter(!is.na(Penalty)) %>%
        mutate(MinorAllele = as.numeric(MinorAllele), MajorAllele = as.numeric(MajorAllele))

    return (df)
}

purityPenaltyPlot <- function(df) {
    ggplot(df) +
        geom_tile(aes(x = MajorAllele, y = MinorAllele,  color = Penalty)) +
        scale_colour_gradientn(colours=c("blue","green","yellow","orange", "red"), limits = c(0, 15)) +
        theme( panel.border = element_blank()) +
        xlab("Major Allele Ploidy") + ylab("Minor Allele Ploidy")
}

ploidy = seq(-1, 8, by = 0.01)
lowPurityDF = purityDataFrame(purityMatrix(0.3, ploidy), ploidy)
highPurityDF = purityDataFrame(purityMatrix(0.7, ploidy), ploidy)

p1 = purityPenaltyPlot(lowPurityDF) + ggtitle("Ploidy Penalty @ 30% Purity") + theme(legend.position = "none")
p2 = purityPenaltyPlot(highPurityDF) + ggtitle("Ploidy Penalty @ 70% Purity")


pfinal = plot_grid(p1,p2, ncol = 2, rel_widths = c(0.9, 1), labels = "AUTO")
save_plot("~/hmf/RPlot/FittedPurityPenalty.png", pfinal, base_width = 8, base_height = 4)



