# Parse the arguments
args <- commandArgs(trailing=T)
bafFile <- args[1]
pcfFile   <- args[2]

library(copynumber)
baf <- read.table(bafFile, header=TRUE)
chromosomeLevels = levels(baf$Chromosome)
chromosomePrefix = ""
if (any(grepl("chr", chromosomeLevels, ignore.case = T))) {
    chromosomePrefix = substr(chromosomeLevels[1], 1, 3)
}

baf <- baf[,c("Chromosome","Position","TumorModifiedBAF")]
baf$Chromosome <- gsub(chromosomePrefix, "", baf$Chromosome, ignore.case = T)
baf.seg<-pcf(baf,verbose=FALSE,gamma=100,kmin=1)
baf.seg$chrom = paste0(chromosomePrefix, baf.seg$chrom)
write.table(baf.seg, file = pcfFile, row.names = F, sep = "\t", quote = F)