# Parse the arguments
args <- commandArgs(trailing=T)
bafFile <- args[1]
pcfFile   <- args[2]
kmin <- 1

library(dplyr)
library(copynumber)
baf <- read.table(bafFile, header=TRUE, stringsAsFactors = T)
chromosomeLevels = levels(baf$chromosome)
chromosomePrefix = ""
if (any(grepl("chr", chromosomeLevels, ignore.case = T))) {
    chromosomePrefix = substr(chromosomeLevels[1], 1, 3)
}

baf <- baf[,c("chromosome","position","tumorModifiedBAF")]
baf$chromosome <- gsub(chromosomePrefix, "", baf$chromosome, ignore.case = T)
baf.seg<-pcf(baf, verbose=FALSE, gamma=100, kmin=kmin)

# copynumber pcf seems to have a bug that causes issue when n.probes == kmin
# we correct it by setting mean to tumorModifiedBAF
baf.seg = left_join(baf.seg, baf, by=c("chrom" = "chromosome", "start.pos" = "position"))
baf.seg$mean = ifelse(baf.seg$n.probes==1, baf.seg$tumorModifiedBAF, baf.seg$mean)

baf.seg = subset(baf.seg, select = -tumorModifiedBAF)
baf.seg$chrom = paste0(chromosomePrefix, baf.seg$chrom)
write.table(baf.seg, file = pcfFile, row.names = F, sep = "\t", quote = F)