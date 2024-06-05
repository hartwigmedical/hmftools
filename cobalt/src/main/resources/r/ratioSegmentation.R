# Parse the arguments
args = commandArgs(trailing = T)
ratioFile = args[1]
column = args[2]
pcfFile = args[3]
gamma = as.numeric(args[4])
kmin = 1

library(dplyr)
ratio = read.table(ratioFile, header = TRUE, stringsAsFactors = T)

chromosomeLevels = levels(ratio$chromosome)
chromosomePrefix = ""
if (any(grepl("chr", chromosomeLevels, ignore.case = T))) {
    chromosomePrefix = substr(chromosomeLevels[1], 1, 3)
}

ratio$Ratio = ratio[, column]
ratio = ratio[ratio$Ratio >= 0,]

if(nrow(ratio) > 0)
{
    library(copynumber)

    ratio$Ratio[ratio$Ratio < 0.001] = 0.001
    ratio$S1 = log2(ratio$Ratio)
    ratio = ratio[! is.nan(ratio$S1),]
    ratio = ratio[, c("chromosome", "position", "S1")]

    ratio$chromosome = gsub(chromosomePrefix, "", ratio$chromosome, ignore.case = T)
    ratio.seg = pcf(ratio, verbose = FALSE, gamma = gamma, kmin = kmin)

    # copynumber pcf seems to have a bug that causes issue when n.probes == 1
    # we correct it by setting mean to tumorModifiedBAF
    ratio.seg = left_join(ratio.seg, ratio, by=c("chrom" = "chromosome", "start.pos" = "position"))
    ratio.seg$mean = ifelse(ratio.seg$n.probes==1, ratio.seg$S1, ratio.seg$mean)

    ratio.seg = subset(ratio.seg, select = -S1)
    ratio.seg$chrom = paste0(chromosomePrefix, ratio.seg$chrom)
    write.table(ratio.seg, file = pcfFile, row.names = F, sep = "\t", quote = F)
} else {
      # write an empty PCF file
      emptyPcf = data.frame(matrix(ncol=7,nrow=0))
      colnames(emptyPcf) = c('sampleID','chrom','arm','start.pos','end.pos','n.probes','mean')
      write.table(emptyPcf, pcfFile, row.names = F, sep = "\t", quote = F)
}
