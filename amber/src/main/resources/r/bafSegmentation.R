# Parse the arguments
args = commandArgs(trailing=T)
bafFile = args[1]
pcfFile = args[2]
kmin = 1

baf = read.table(bafFile, header=TRUE, stringsAsFactors = T)

library(dplyr)

if(nrow(baf) > 0)
{
    library(copynumber)

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
    write.table(baf.seg, pcfFile, row.names = F, sep = "\t", quote = F)
} else {
    # write an empty PCF file
    emptyPcf = data.frame(matrix(ncol=7,nrow=0))
    colnames(emptyPcf) = c('sampleID','chrom','arm','start.pos','end.pos','n.probes','mean')
    write.table(emptyPcf, pcfFile, row.names = F, sep = "\t", quote = F)
}

