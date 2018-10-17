# Parse the arguments
args <- commandArgs(trailing=T)
bafFile <- args[1]
pcfFile   <- args[2]
pcfFile1  <- paste(args[2], "1", sep="")

library(copynumber)
baf <- read.table(bafFile, header=TRUE)
baf <- baf[,c("Chromosome","Position","TumorModifiedBAF")]
baf$Chromosome <- gsub("chr", "", baf$Chromosome, ignore.case = T)
baf.seg<-pcf(baf,verbose=FALSE,gamma=100,kmin=1,save.res = TRUE, file.names = c(pcfFile1, pcfFile))