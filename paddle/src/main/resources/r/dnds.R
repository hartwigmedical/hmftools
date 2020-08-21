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
refdb="/Users/jon/hmf/repos/hmftools/paddle/src/main/resources/r/HmfRefCDS.RData"
load(refdb)


files = list.files("/Users/jon/hmf/analysis/dnds/somatics/", full.names = T, pattern = "exonic.somatics.tsv")[1:1000]
allSomatics = data.frame(stringsAsFactors = F)
for (file in files) {
  somatics = read.table(file = file, sep = "\t", header = T, stringsAsFactors = F)
  exonicSomatics = exonic_somatics(somatics, gr_genes)
  allSomatics = bind_rows(allSomatics, exonicSomatics)
}

rm(somatics)
rm(exonicSomatics)

#save(allSomatics, file = "/Users/jon/hmf/analysis/dnds/somatics/allSomatics.RData")


newgenes = sapply(RefCDS, function(x) x$gene_name) # New list of genes
oldgenes = sapply(strsplit(newgenes,split=":"), function(x) x[1])
cv = covs[oldgenes,]
rownames(cv) = newgenes
kc = newgenes[oldgenes %in% known_cancergenes]

output = dndscv(allSomatics, refdb=refdb, kc=kc, cv=cv, stop_loss_is_nonsense = TRUE)
HmfRefCDSCv <- output$sel_cv
write.table(file = "/Users/jon/hmf/repos/hmftools/paddle/src/main/resources/HmfRefCDSCv.new.tsv", HmfRefCDSCv, quote = F, row.names = F, sep = "\t")



dndsFilteredAnnotatedMutations = output$annotmuts

head(allSomatics)

