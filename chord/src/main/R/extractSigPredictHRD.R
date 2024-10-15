#!/usr/bin/env Rscript

options(stringsAsFactors=F) # to avoid invalid factor level warning

args <- commandArgs(TRUE)

outDir <- args[1]
sampleName <- args[2]
snvIndVcf <- args[3]
svVcf <- args[4]
refGenomeVsn <- args[5] # V37 or V38
chordToolDir <- args[6]

sigOutTxt <- paste0(outDir, '/', sampleName, '_chord_signatures.txt')
prdOutTxt <- paste0(outDir, '/', sampleName, '_chord_prediction.txt')

LOG_LEVEL <- list(
  FATAL="FATAL",
  ERROR="ERROR",
  WARN ="WARN ",
  INFO ="INFO ",
  DEBUG="DEBUG",
  TRACE="TRACE"
)

logMessage <- function(..., log.level = LOG_LEVEL$INFO){
  current_time <- format(Sys.time(), "%H:%H:%OS3")
  message(current_time, " [R] [", log.level, "] ", ...)
}

logMessage("START CHORD signature extraction and HRD prediction")
setwd(outDir)

logMessage("CHORD Settings: ")
logMessage("  Sample name: ", sampleName)
logMessage("  Somatic SNV/IND vcf: ", snvIndVcf)
logMessage("  Somatic SV vcf: ", svVcf)
logMessage("  Ref Genome Version: ", refGenomeVsn)
logMessage("  Output dir: ", outDir)
logMessage("  Signature out file: ", sigOutTxt)
logMessage("  Prediction out file: ", prdOutTxt)
logMessage("  CHORD tool dir: ", chordToolDir)

## Load packages
suppressPackageStartupMessages(library('randomForest'))

if(!is.na(chordToolDir)){
  logMessage("Loading mutSigExtractor and CHORD from source dir: ", chordToolDir)
  suppressPackageStartupMessages(library('devtools'))
  suppressPackageStartupMessages(devtools::load_all(paste0(chordToolDir, '/mutSigExtractor')))
  suppressPackageStartupMessages(devtools::load_all(paste0(chordToolDir, '/CHORD')))
} else {
  logMessage("Loading mutSigExtractor and CHORD from library: ", paste(.libPaths(), collapse=", "))
  ## CHORD and mutSigExtractor can be installed like so:
  ## Rscript -e "install.packages('hmftools/chord/src/main/R/mutSigExtractor', repos = NULL, type='source')"
  suppressPackageStartupMessages(library('mutSigExtractor'))
  suppressPackageStartupMessages(library('CHORD'))
}

## Convert genome name to BSGenome name
if (refGenomeVsn == "V37") {
  suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
  refGenome <- BSgenome.Hsapiens.UCSC.hg19
} else if (refGenomeVsn == "V38") {
  suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  refGenome <- BSgenome.Hsapiens.UCSC.hg38
} else {
  logMessage("Unsupported ref genome version: ", refGenomeVsn," (should be V37 or V38)", log.level=LOG_LEVEL$ERROR)
  stop()
}

## Extract signature counts and predict HRD
logMessage("Performing CHORD signature extraction")
signatures <- CHORD::extractSigsChord(
  vcf.snv = snvIndVcf,
  vcf.indel = snvIndVcf,
  vcf.sv = svVcf,
  sample.name = sampleName,
  sv.caller = "esvee",
  vcf.filters=list(snv="PASS", indel="PASS", sv="PASS"),
  ref.genome=refGenome
)
logMessage("Performing CHORD HRD prediction")
prediction <- CHORD::chordPredict(
  signatures,
  rf.model = CHORD,
  hrd.cutoff = 0.5
)

## Output
logMessage("Writing output file: ", sigOutTxt)
write.table(signatures, file=sigOutTxt, sep="\t")

logMessage("Writing output file: ", prdOutTxt)
write.table(prediction, file=prdOutTxt, sep="\t", quote=FALSE, row.names=FALSE)

logMessage("FINISHED CHORD signature extraction and HRD prediction")
