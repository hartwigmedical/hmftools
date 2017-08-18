# Parse the arguments
args <- commandArgs(trailing=T)
ratioFile <- args[1]
pcfFile   <- args[2]
pcfFile1  <- paste(args[2], "1", sep="")

library(copynumber)
ratio <- read.table(ratioFile, header=TRUE)
ratio <- ratio[ratio$Ratio>=0,]
ratio$Ratio[ratio$Ratio<0.001] <- 0.001
ratio$S1 = log2(ratio$Ratio)
ratio <- ratio[!is.nan(ratio$S1),]
ratio <- ratio[,c("Chromosome","Position","S1")]
ratio.seg<-pcf(ratio,verbose=FALSE,gamma=100,kmin=1,save.res = TRUE, file.names = c(pcfFile1, pcfFile))