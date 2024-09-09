#!/usr/bin/env Rscript
## Run on hpc
library(devtools)
load_all('/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/hmfGeneAnnotation/')

args <- commandArgs(trailingOnly=T)

detGeneStatuses(
   out.dir = args[1],
   input.file.paths = c(germ_vcf=args[2], som_vcf=args[3], gene_cnv=args[4], cnv=args[5]),
   sample.name = args[6],
   java.path = args[7],
   bed.file = args[8]
)

