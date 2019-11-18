# Parse the arguments
args <- commandArgs(trailing = T)
ratioFile <- args[1]
column <- args[2]
pcfFile <- args[3]

library(copynumber)
ratio <- read.table(ratioFile, header = TRUE, stringsAsFactors = T)

chromosomeLevels = levels(ratio$chromosome)
chromosomePrefix = ""
if (any(grepl("chr", chromosomeLevels, ignore.case = T))) {
    chromosomePrefix = substr(chromosomeLevels[1], 1, 3)
}

ratio$Ratio = ratio[, column]
ratio <- ratio[ratio$Ratio >= 0,]
ratio$Ratio[ratio$Ratio < 0.001] <- 0.001
ratio$S1 = log2(ratio$Ratio)
ratio <- ratio[! is.nan(ratio$S1),]
ratio <- ratio[, c("chromosome", "position", "S1")]

ratio$chromosome <- gsub(chromosomePrefix, "", ratio$chromosome, ignore.case = T)
ratio.seg <- pcf(ratio, verbose = FALSE, gamma = 100, kmin = 1)
ratio.seg$chrom = paste0(chromosomePrefix, ratio.seg$chrom)
write.table(ratio.seg, file = pcfFile, row.names = F, sep = "\t", quote = F)
