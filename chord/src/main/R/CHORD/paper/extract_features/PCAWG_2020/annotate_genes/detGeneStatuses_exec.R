#!/usr/bin/env Rscript
## Run on hpc
library(devtools)
load_all('/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/')

args <- commandArgs(trailingOnly=T)

detGeneStatuses(
   out.dir = args[1],
   input.file.paths = c(germ_vcf=args[2], som_vcf=args[3], gene_cnv=NA, cnv=args[4]),
   sample.name = args[5],
   java.path = args[6],
   bed.file = args[7],
   exons.bed.file = args[8],
   chrom.arm.split.method='universal',
   do.filter.vcf=T, do.snpeff.ann=T, do.det.reading.frame=F,
   ignore.chroms='Y'
)

