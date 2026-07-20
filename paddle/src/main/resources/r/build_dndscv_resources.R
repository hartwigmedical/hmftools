ARGS <- optparse::parse_args(optparse::OptionParser(option_list = list(
      optparse::make_option("--ensembl_data_dir", type = "character", help = "Dir containing ensembl_trans_exon_data.csv / ensembl_gene_data.csv"),
      optparse::make_option("--driver_gene_panel", type = "character", help = "Path to DriverGenePanel TSV"),
      optparse::make_option("--genome_fasta", type = "character", help = "Path to reference genome FASTA"),
      optparse::make_option("--dndscv_package_dir", type = "character", help = "Path to dndscv package checkout (for bundled covariates data)"),
      optparse::make_option("--output_dir", type = "character", help = "Dir to write output resources to")
   )
))

ENSEMBL_DATA_DIR       <- ARGS$ensembl_data_dir
DRIVER_GENE_PANEL_PATH <- ARGS$driver_gene_panel
GENOME_FASTA_PATH      <- ARGS$genome_fasta
DNDSCV_PACKAGE_DIR     <- ARGS$dndscv_package_dir
OUTPUT_DIR             <- ARGS$output_dir

if(FALSE){
   ENSEMBL_DATA_DIR <- "/Users/lnguyen/Hartwig/repos/pipeline-resources/hmftools/common/ensembl_data/"
   DRIVER_GENE_PANEL_PATH <- "/Users/lnguyen/Hartwig/repos/pipeline-resources/hmftools/common/DriverGenePanel.38.tsv"
   GENOME_FASTA_PATH <- "/Users/lnguyen/Hartwig/repos/pipeline-resources/reference_genome/38/GRCh38_masked_exclusions_alts_hlas.fasta"
   DNDSCV_PACKAGE_DIR <- "/Users/lnguyen/Hartwig/repos/dndscv/"
   OUTPUT_DIR <- "/Users/lnguyen/Hartwig/experiments/paddle/20260714_run_dndscv/resources/"
}

## Coding exon coords ================================

## CDKN2A encodes two unrelated proteins (p16INK4a, p14ARF) from overlapping
## reading frames at the same locus. dndscv's own hg19/hg38 RefCDS
## (and its bundled covariates) model these as two separate "genes",
## CDKN2A.p16INK4a and CDKN2A.p14arf, rather than one CDKN2A gene, so both
## transcripts are kept here instead of just the canonical one.
CDKN2A_TRANSCRIPTS <- c(ENST00000304494 = "CDKN2A.p16INK4a", ENST00000579755 = "CDKN2A.p14arf")

get_exon_coords <- function(){
   
   message("Getting exon coordinates")
   
   exons <- read.csv(paste0(ENSEMBL_DATA_DIR,"/ensembl_trans_exon_data.csv"), stringsAsFactors = FALSE)
   genes <- read.csv(paste0(ENSEMBL_DATA_DIR,"/ensembl_gene_data.csv"), stringsAsFactors = FALSE)
   
   ## Handle CDKN2A
   exon_coords <- exons |>
      
      ## Canonical transcripts (+ the CDKN2A isoforms), with gene name/chromosome attached
      dplyr::filter(
         BioType == "protein_coding",
         TransId == CanonicalTranscriptId | TransName %in% names(CDKN2A_TRANSCRIPTS)
      ) |>
      dplyr::left_join(dplyr::select(genes, GeneId, GeneName, Chromosome), by = "GeneId") |>
      
      ## Split CDKN2A into its two isoform-specific gene names
      dplyr::mutate(
         GeneName = ifelse(TransName %in% names(CDKN2A_TRANSCRIPTS), CDKN2A_TRANSCRIPTS[TransName], GeneName)
      )
   
   ## Clip exon boundaries to the CDS region
   exon_coords <- exon_coords |>
      dplyr::mutate(
         CodingStart = as.integer(CodingStart),
         CodingEnd   = as.integer(CodingEnd),
         chr_coding_start = pmax(ExonStart, CodingStart),
         chr_coding_end   = pmin(ExonEnd, CodingEnd)
      )
   
   ## Drop exons that end up with no coding sequence
   exon_coords <- exon_coords |>
      dplyr::filter(!is.na(chr_coding_start), chr_coding_start <= chr_coding_end) |>
      dplyr::mutate(coding_length = chr_coding_end - chr_coding_start + 1L)
   
   ## Sort exons into transcription order (plus strand ascending, minus strand descending)
   exon_coords <- exon_coords |>
      dplyr::arrange(TransId, ifelse(Strand == 1L, ExonStart, -ExonStart))
   
   ## Compute each transcript's CDS-relative coordinates
   exon_coords <- exon_coords |>
      dplyr::group_by(TransId) |>
      dplyr::mutate(
         cds_start  = cumsum(coding_length) - coding_length + 1L,
         cds_end    = cumsum(coding_length),
         cds_length = sum(coding_length)
      ) |>
      dplyr::ungroup()
   
   ## Select/rename the final columns and sort genome-wide by chromosome and position
   exon_coords <- exon_coords |>
      dplyr::transmute(
         gene.id          = GeneId,
         gene.name        = GeneName,
         cds.id           = TransName,
         chr              = Chromosome,
         chr.coding.start = chr_coding_start,
         chr.coding.end   = chr_coding_end,
         cds.start        = cds_start,
         cds.end          = cds_end,
         length           = cds_length,
         strand           = Strand,
         exon.rank        = ExonRank
      ) |>
      dplyr::arrange(factor(chr, c(as.character(1:22), "X", "Y")), chr.coding.start)
   
   return(exon_coords)
}
EXON_COORDS <- get_exon_coords()


## Build ref CDS, used by dndscv ================================
build_dndscv_ref_cds_data <- function(){
   
   message("Running dndscv::buildref")
   
   ref_cds_rda <- paste0(OUTPUT_DIR, "/RefCDS_GRCh38.rda")
   
   if(!file.exists(ref_cds_rda)){
      EXON_COORDS_tsv <- paste0(OUTPUT_DIR, "/exon_coords.GRCh38.tsv")
      write.table(EXON_COORDS, EXON_COORDS_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
      dndscv::buildref(cdsfile = EXON_COORDS_tsv, genomefile = GENOME_FASTA_PATH, outfile = ref_cds_rda, excludechrs = "MT")
   }
}
build_dndscv_ref_cds_data()


## Make bed file to slice VCFs ================================
create_exon_coords_bed <- function(){
   
   message("Creating exons bed file")
   
   load(paste0(OUTPUT_DIR, "/RefCDS_GRCh38.rda"))  # RefCDS, gr_genes
   
   ## Flatten contiguous rows
   exon_coords_gr <- GenomicRanges::reduce(gr_genes, with.revmap = TRUE)                                                                                                                   
   
   ## Restore gene names
   revmap <- S4Vectors::mcols(exon_coords_gr)$revmap
   gene_names <- IRanges::extractList(gr_genes$names, revmap)
   
   gene_names <- sapply(gene_names, function(genes){
      genes <- unique(genes)
      
      if(length(genes)>1){
         paste(genes, collapse=';')
      } else {
         genes
      }
   })
   exon_coords_gr$gene <- gene_names
   
   ## Sort by coordinates
   exon_coords_gr <- exon_coords_gr |> GenomeInfoDb::sortSeqlevels() |> BiocGenerics::sort()
   
   ## Write bed file
   exon_coords_bed <- exon_coords_gr |> 
      as.data.frame() |>
      dplyr::transmute(
         chr   = seqnames,
         start = start - 1L, # convert 1-based -> 0-based
         end   = end,
         gene  = gene
      )
   
   write.table(
      exon_coords_bed, paste0(OUTPUT_DIR, "/exon_coords.GRCh38.bed"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
   )
}
create_exon_coords_bed()


## Fix covariates gene names ================================

create_dndscv_covariates_matrix <- function(){
   
   message("Creating dndscv covariates matrix")
   
   ## Our RefCDS uses current HGNC gene symbols, but the covariates matrix (used
   ## by dndscv for its negative-binomial background-rate regression) was built
   ## on an older HGNC annotation, so some RefCDS genes are missing under their
   ## current symbol. For each missing gene, look up its previous/alias symbols
   ## in HGNC's complete gene set and, if one of them is already a row in the
   ## covariates matrix, rename that row to the current symbol. CDKN2A is
   ## excluded: dndscv models it as two isoform-specific genes (CDKN2A.p16INK4a /
   ## CDKN2A.p14arf) by design, not as a naming-vintage mismatch, and both names
   ## already exist in the covariates matrix (see buildref of RefCDS above).
   
   hgnc_file <- paste0(OUTPUT_DIR, "/hgnc_complete_set.txt")
   if (!file.exists(hgnc_file)) {
      download.file("https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt", hgnc_file)
   }
   hgnc <- read.delim(hgnc_file, quote = "\"", stringsAsFactors = FALSE, na.strings = "")
   
   load(paste0(DNDSCV_PACKAGE_DIR,"/data/covariates_hg19_hg38_epigenome_pcawg.rda"))  # covs
   
   refcds_genes  <- unique(EXON_COORDS$gene.name)
   missing_genes <- setdiff(refcds_genes, rownames(covs))
   missing_genes <- setdiff(missing_genes, names(CDKN2A_TRANSCRIPTS))
   
   ## For a gene missing from the covariates matrix, find a previous/alias
   ## symbol (per HGNC) that already exists as a row in the covariates matrix
   find_covariate_alias <- function(gene, hgnc, covariate_names) {
      hgnc_row <- hgnc[hgnc$symbol == gene, ]
      if (nrow(hgnc_row) == 0) return(NA_character_)
      
      candidates <- unlist(strsplit(c(hgnc_row$prev_symbol[1], hgnc_row$alias_symbol[1]), "|", fixed = TRUE))
      candidates <- candidates[!is.na(candidates)]
      
      match <- candidates[candidates %in% covariate_names]
      if (length(match) == 0) return(NA_character_)
      match[1]
   }
   
   covariate_aliases <- lapply(missing_genes, find_covariate_alias, hgnc = hgnc, covariate_names = rownames(covs))
   names(covariate_aliases) <- missing_genes
   covariate_aliases <- unlist(covariate_aliases[!is.na(covariate_aliases)])
   
   cat(sprintf(
      "%d of %d genes missing from the covariates matrix resolved via HGNC prior/alias symbols; %d remain unresolved:\n",
      length(covariate_aliases), length(missing_genes), length(missing_genes) - length(covariate_aliases)
   ))
   print(setdiff(missing_genes, names(covariate_aliases)))
   
   ## Rename the resolved covariates rows to the current gene symbol
   rownames(covs)[match(covariate_aliases, rownames(covs))] <- names(covariate_aliases)
   
   write.table(
      covs,
      gzfile(paste0(OUTPUT_DIR, "/covariates.GRCh38.tsv.gz")),
      sep = '\t', quote = F, row.names = T
   )
}
create_dndscv_covariates_matrix()


## Prepare known cancer genes (kc) ================================
create_dndscv_genes_list <- function(){
   ## dndscv's `kc` argument is the list of a-priori known cancer genes to
   ## exclude from the background/indel mutation-rate model. Use Hartwig's own
   ## curated driver gene panel (already used by Purple/Orange) instead of
   ## dndscv's generic bundled Cancer Gene Census snapshot, for consistency with
   ## the rest of the driver-calling pipeline.
   driver_gene_panel <- read.delim(DRIVER_GENE_PANEL_PATH)
   driver_genes <- driver_gene_panel$gene
   
   ## The panel lists a single CDKN2A gene, but our RefCDS splits it into two
   ## isoform-specific genes (CDKN2A.p16INK4a / CDKN2A.p14arf); substitute so
   ## both are excluded from the background model. This mirrors what dndscv does
   ## internally for its own bundled hg19/hg38 refdb, which isn't triggered here
   ## since we use a custom refdb.
   if ("CDKN2A" %in% driver_genes) {
      driver_genes <- unique(c(setdiff(driver_genes, "CDKN2A"), "CDKN2A.p16INK4a", "CDKN2A.p14arf"))
   }
   
   writeLines(driver_genes, paste0(OUTPUT_DIR, "/driver_genes.GRCh38.txt"))
}
create_dndscv_genes_list()


