
args <- commandArgs(trailing = T)

work_dir <- args[1]

somatics_input_dir <- paste0(work_dir, "/somatics")
mutation_load_input_file <- paste0(work_dir, "/mutationalLoad.tsv")

all_somatics_output_file <- paste0(work_dir, "/allSomatics.RData")
ref_cdscv_output_rdata <- paste0(work_dir, "/HmfRefCDSCv.RData")
ref_cdscv_output_tsv <- paste0(work_dir, "/HmfRefCDSCv.tsv")
dnds_mutations_output_tsv <- paste0(work_dir, "/DndsMutations.tsv")

library(GenomicRanges)
library(dplyr)
library(tidyr)
library(dndscv)

exonic_somatics <- function(somatics, gr_genes) {
  gr_muts = GRanges(somatics$chromosome, IRanges(somatics$position,somatics$position + nchar(somatics$ref) - 1))
  ol = as.matrix(findOverlaps(gr_muts, gr_genes, type="any", select="all"))
  return (somatics[unique(ol[, 1]), ])
}

data("cancergenes_cgc81", package="dndscv")
data("covariates_hg19", package="dndscv")
refdb="/data/resources/custom/dnds/HmfRefCDS.RData"
load(refdb)

summary = read.table(file = mutation_load_input_file, header = T)

files = list.files(somatics_input_dir, full.names = T, pattern = "exonic.somatics.tsv")
allSomatics = data.frame(stringsAsFactors = F)
for (i in 1:length(files)) {
  file = files[i]
  cat("Processing", file, "\n")

  somatics = read.table(file = file, sep = "\t", header = T, stringsAsFactors = F)
  exonicSomatics = exonic_somatics(somatics, gr_genes)
  exonicSomatics$chromosome = as.character(exonicSomatics$chromosome)
  allSomatics = bind_rows(allSomatics, exonicSomatics)
}

rm(somatics)
rm(exonicSomatics)

save(allSomatics, file = all_somatics_output_file)

############### Dnds Results
newgenes = sapply(RefCDS, function(x) x$gene_name) # New list of genes
oldgenes = sapply(strsplit(newgenes,split=":"), function(x) x[1])
cv = covs[oldgenes,]
rownames(cv) = newgenes
kc = newgenes[oldgenes %in% known_cancergenes]

output = dndscv(allSomatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)
HmfRefCDSCv <- output$sel_cv
save(HmfRefCDSCv, file = ref_cdscv_output_rdata)
write.table(HmfRefCDSCv, file = ref_cdscv_output_tsv, quote = F, row.names = F, sep = "\t")

############### Dnds Mutations
unfilteredSomatics = allSomatics %>% select(sample, chr = chromosome, pos = position, ref, alt, worstCodingEffect, canonicalCodingEffect, repeatCount, biallelic, hotspot)
unfilteredOutput = dndscv(unfilteredSomatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE, max_muts_per_gene_per_sample = 30000000, max_coding_muts_per_sample = 30000000)
dndsUnfilteredAnnotatedMutations = unfilteredOutput$annotmuts

names(dndsUnfilteredAnnotatedMutations) <- c("sample", "chr","pos", "ref", "mut", "worstCodingEffect", "canonicalCodingEffect", "repeatCount", "biallelic", "hotspot", "geneind", "gene", "ref_cod", "mut_cod", "strand", "ref3_cod", "mut3_cod", "aachange", "ntchange","impact",  "pid")
dndsUnfilteredAnnotatedMutations = dndsUnfilteredAnnotatedMutations %>% select(sample, chr, pos, ref, mut, worstCodingEffect, canonicalCodingEffect, repeatCount, biallelic, hotspot, gene, impact)
write.table(dndsUnfilteredAnnotatedMutations, file = dnds_mutations_output_tsv, quote = F, row.names = F, sep = "\t")
