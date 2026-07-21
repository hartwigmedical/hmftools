ARGS <- optparse::parse_args(optparse::OptionParser(option_list = list(
      optparse::make_option("--cohort_mutations", type = "character", help = "Path to cohort mutations TSV (SampleId, Chromosome, Position, Ref, Alt)"),
      optparse::make_option("--dndscv_resources_dir", type = "character", help = "Dir containing RefCDS_GRCh38.rda / covariates.GRCh38.tsv.gz / driver_genes.GRCh38.txt"),
      optparse::make_option("--dndscv_package_dir", type = "character", help = "Path to dndscv package repo"),
      optparse::make_option("--output_dir", type = "character", help = "Dir to write dndscv output to")
   )
))

COHORT_MUTATIONS_PATH <- ARGS$cohort_mutations
DNDSCV_RESOURCES_DIR  <- ARGS$dndscv_resources_dir
DNDSCV_PACKAGE_DIR    <- ARGS$dndscv_package_dir
OUTPUT_DIR            <- ARGS$output_dir

if(FALSE){
   COHORT_MUTATIONS_PATH <- "/Users/lnguyen/Hartwig/experiments/paddle/20260714_run_dndscv/sample_data/dnds_cohort_variants.tsv.gz"
   OUTPUT_DIR <- "/Users/lnguyen/Hartwig/experiments/paddle/20260714_run_dndscv/dndscv_output/"
   DNDSCV_RESOURCES_DIR <- "/Users/lnguyen/Hartwig/experiments/paddle/20260714_run_dndscv/resources/"
   DNDSCV_PACKAGE_DIR <- "/Users/lnguyen/Hartwig/repos/dndscv/"
}


devtools::load_all(DNDSCV_PACKAGE_DIR)

## Load resources
load(paste0(DNDSCV_RESOURCES_DIR,"/RefCDS_GRCh38.rda")) # RefCDS, gr_genes
covariates <- read.delim(paste0(DNDSCV_RESOURCES_DIR,"/covariates.GRCh38.tsv.gz"))
driver_genes <- read.delim(paste0(DNDSCV_RESOURCES_DIR,"/driver_genes.GRCh38.txt"), header = FALSE)[[1]]

## Load cohort mutations
mutations <- read.delim(COHORT_MUTATIONS_PATH)
mutations <- mutations[, c("SampleId", "Chromosome", "Position", "Ref", "Alt")]

## Run dndscv
dndscv_out <- dndscv::dndscv(mutations=mutations, refdb=RefCDS, kc=driver_genes, cv=covariates)

dir.create(OUTPUT_DIR)

write.table(
   dndscv_out$sel_cv, 
   gzfile(paste0(OUTPUT_DIR,"/dndscv_output.sel_cv.tsv.gz")),
   sep = '\t', row.names = F, quote = F
)

saveRDS(dndscv_out, paste0(OUTPUT_DIR,"/dndscv_output.rds"))
