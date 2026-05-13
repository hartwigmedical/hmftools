library(stringr)
library(dplyr)
library(GenomicRanges)

library(ggplot2)
library(patchwork)

theme_set(
   theme_bw() +
   theme(
      panel.grid.minor = element_blank()
   )
)

## =============================
## Config
## =============================

#' Usage example:
#'
#' Rscript compare_cohorts.R \
#' --reference_input_dir /Users/lnguyen/Hartwig/experiments/wigits_qc/analysis/20260414_panel_metrics/hmf_data \
#' --target_input_dir /Users/lnguyen/Hartwig/experiments/wigits_qc/analysis/20260414_panel_metrics/panel_data/msk \
#' --panel_bed_file /Users/lnguyen/Hartwig/resources/common-resources-public/panel/msk/panel_definition.msk.37.bed.gz \
#' --driver_gene_panel /Users/lnguyen/Hartwig/resources/common-resources-public/panel/msk/driver_genes.msk.37.tsv \
#' --output_dir /Users/lnguyen/Hartwig/experiments/wigits_qc/analysis/20260414_panel_metrics/output

args <- local({
   parser <- argparse::ArgumentParser()
   
   parser$add_argument("--reference_input_dir", help="Directory containing reference cohort data")
   parser$add_argument("--target_input_dir", help="Directory containing the target cohort data")
   
   parser$add_argument("--reference_name", help="Reference cohort name. Used in plots and output file prefixes", default="reference")
   parser$add_argument("--target_name", help="Target cohort name. Used in plots and output file prefixes", default="target")
   parser$add_argument("--output_dir", help="Output directory")
   
   parser$add_argument("--reference_tables_dir", help="Directory containing reference cohort merged TSVs. Defaults: <output_dir>/tables")
   parser$add_argument("--target_tables_dir", help="Directory containing target cohort merged TSVs. Defaults: <output_dir>/tables")
   
   parser$add_argument("--panel_bed_file", help="Path to the panel bed file")
   parser$add_argument("--driver_gene_panel", help="Path to the driver gene panel file")
   
   parser$add_argument("--log_level", help="Log level", default = "DEBUG")
   
   parser$parse_args()
})

if(is.null(args$reference_tables_dir)) args$reference_tables_dir <- file.path(args$output_dir, "tables")
if(is.null(args$target_tables_dir)) args$target_tables_dir <- file.path(args$output_dir, "tables")

## =============================
## Logging / common functions
## =============================

LOG_LEVEL <- list(
   TRACE = list(name = "TRACE", severity = 1),
   DEBUG = list(name = "DEBUG", severity = 2),
   INFO  = list(name = "INFO" , severity = 3),
   WARN  = list(name = "WARN" , severity = 4),
   ERROR = list(name = "ERROR", severity = 5),
   FATAL = list(name = "FATAL", severity = 6)
)

log_message <- function(log_level, fmt, ...){
   
   current_time <- format(Sys.time(), "%H:%H:%OS3")
   
   log_message <- sprintf("%s [R] [%-5s] %s", current_time, log_level$name, sprintf(fmt, ...))
   
   if(log_level$severity >= LOG_LEVEL[[LOG_LEVEL$ERROR$name]]$severity)
      stop(log_message)
   
   if(log_level$severity >= LOG_LEVEL[[toupper(args$log_level)]]$severity)
      message(log_message)
}

LOGGER <- list(
   trace = function(fmt, ...){ log_message(LOG_LEVEL$TRACE, fmt, ...) },
   debug = function(fmt, ...){ log_message(LOG_LEVEL$DEBUG, fmt, ...) },
   info  = function(fmt, ...){ log_message(LOG_LEVEL$INFO , fmt, ...) },
   warn  = function(fmt, ...){ log_message(LOG_LEVEL$WARN , fmt, ...) },
   error = function(fmt, ...){ log_message(LOG_LEVEL$ERROR, fmt, ...) },
   fatal = function(fmt, ...){ log_message(LOG_LEVEL$FATAL, fmt, ...) }
)

## =============================
## Resources
## =============================

PANEL_BED <- local({
   
   LOGGER$info("Loading panel bed file: %s", args$panel_bed_file)
   
   maybe_header <- read.delim(args$panel_bed_file, nrows = 1, header = FALSE)
   has_header <- grepl("\\w+", maybe_header[,2])
   bed <- read.delim(args$panel_bed_file, header = has_header)
   
   bed <- bed[1:3]
   colnames(bed) <- c("Chromosome", "PosStart", "PosEnd")
   
   return(bed)
})

DRIVER_GENES <- local({
   
   LOGGER$info("Loading driver gene panel: %s", args$driver_gene_panel)
   
   driver_gene_panel <- read.delim(args$driver_gene_panel)
   
   driver_gene_panel <- lapply(driver_gene_panel, function(column){
      if(all(column %in% c("true","false"))){
         column <- as.logical(column)
      }
      return(column)
   }) %>% as.data.frame()
   
   driver_genes <- list()
   
   driver_genes$AMP <- driver_gene_panel %>% dplyr::filter(reportAmplification) %>% pull(gene)
   driver_genes$DEL <- driver_gene_panel %>% dplyr::filter(reportDeletion) %>% pull(gene)
   driver_genes$CNV <- sort(c(driver_genes$AMP, driver_genes$DEL))
   
   driver_genes$RNA <- NULL 
   if(all(c("reportHighExpression", "reportLowExpression") %in% colnames(driver_gene_panel))){
      driver_genes$RNA <- driver_gene_panel %>% dplyr::filter(reportHighExpression | reportLowExpression) %>% pull(gene)   
   }
   
   return(driver_genes)
})

## =============================
## Config
## =============================

COHORT_TYPE <- list(
   
   REFERENCE = list(
      name = args$reference_name, 
      input_dir = args$reference_input_dir,
      tables_dir = args$reference_tables_dir
   ),
   
   TARGET = list(
      name = args$target_name, 
      input_dir = args$target_input_dir,
      tables_dir = args$target_tables_dir
   )
   
)

DEFAULT_COVERAGE_THRESHOLDS <- c(50, 100, 200, 400)

## =============================
## Loader functions 
## =============================

if(!dir.exists(args$reference_tables_dir)) dir.create(args$reference_tables_dir, recursive = TRUE)
if(!dir.exists(args$target_tables_dir)) dir.create(args$target_tables_dir, recursive = TRUE)

merge_tsvs <- function(cohort_type, pattern, output_file_suffix, read_func = NULL, progress_interval = 500){
   
   output_path <- file.path(
      cohort_type$tables_dir, 
      paste0(cohort_type$name, output_file_suffix)
   )
   
   if(!is.null(output_path) && file.exists(output_path)){
      LOGGER$info("Loading existing output file: %s", output_path)
      return(read.delim(output_path))
   }

   directory <- cohort_type$input_dir
   
   paths <- list.files(directory, pattern, full.names = TRUE, recursive = TRUE)
   if(length(paths) == 0){
      LOGGER$warn("No files with pattern(%s) found at: %s", pattern, directory)
      return(NULL)
   }

   LOGGER$info("Loading files with pattern '%s' at: %s", pattern, directory)

   files_info <- data.frame(
      SampleId = basename(paths) %>% strsplit('.', fixed = TRUE) %>% sapply(`[[`, 1),
      Path = paths
   )

   if(is.null(read_func)){
      read_func <- read.delim
   }
   
   df_merged <- lapply(1:nrow(files_info), function(i){
      path <- files_info$Path[i]
      sample_id <- files_info$SampleId[i]

      if(i %% progress_interval == 0){
         LOGGER$debug("  Loaded %s/%s files. Current file: %s", i, nrow(files_info), path)
      }

      df <- read_func(path)
      df <- df %>% dplyr::mutate(SampleId = sample_id, .before = 1)
      return(df)
   })

   df_merged <- do.call(rbind, df_merged)
   df_merged$SampleId <- factor(df_merged$SampleId, files_info$SampleId)

   if(!is.null(output_path)){
      write.table(df_merged, gzfile(output_path), sep = "\t", quote = F, row.names = F)
   }

   return(df_merged)
}

## --------------------------------
## BamMetrics
## --------------------------------

merge_tsvs.bam_metric.summary <- function(cohort_type){ 
   merge_tsvs(cohort_type, pattern="*.bam_metric.summary.tsv(.gz)*$", output_file_suffix=".bam_metric.summary.tsv.gz") 
}

merge_tsvs.bam_metric.frag_length <- function(cohort_type){
   
   read_func <- function(path){
      
      df <- read.delim(path)
      
      which_median <- function(x){
         mid_point <- sum(x) / 2
         min(which(cumsum(x) >= mid_point))
      }
      
      data.frame(MedianFragmentLength = df$FragmentLength[which_median(df$Count)])
   }
   
   merge_tsvs(cohort_type, pattern = "*.bam_metric.frag_length.tsv(.gz)*$", output_file_suffix = ".bam_metric.frag_length.tsv.gz", read_func = read_func)
}

merge_tsvs.bam_metric.coverage <- function(cohort_type, coverage_thresholds=DEFAULT_COVERAGE_THRESHOLDS){
   
   read_func <- function(path){
      df <- read.delim(path)
      
      bases_above_coverage <- sapply(coverage_thresholds, function(coverage_threshold){ sum(df$Count[df$Coverage > coverage_threshold]) })
      total_bases <- sum(df$Count)
      prop_bases_above_coverage <- bases_above_coverage / total_bases
      
      data.frame(CoverageAbove = coverage_thresholds, PropBasesAboveCoverage = prop_bases_above_coverage)
   }
   
   merge_tsvs(cohort_type, pattern="*.bam_metric.coverage.tsv(.gz)*$", output_file_suffix = ".bam_metric.coverage.tsv.gz", read_func=read_func)
}

merge_tsvs.bam_metric.gene_coverage <- function(cohort_type, coverage_thresholds=DEFAULT_COVERAGE_THRESHOLDS){
   
   read_func <- function(path){
      
      df <- read.delim(path)
      
      ## Get depth matrix
      depths <- df %>% dplyr::select(dplyr::starts_with("DR")) %>% as.matrix()
      rownames(depths) <- df$GeneName
      
      ## Get depth bins from header
      bin_min <- stringr::str_extract(colnames(depths), "_(\\d+)", group = 1) %>% as.numeric()
      
      ## Calc proportion of bases in gene above X coverage
      total_gene_bases <- rowSums(depths)
      
      output <- lapply(coverage_thresholds, function(coverage_threshold){
         #coverage_threshold=50
         bins_above_coverage <- depths[,bin_min >= coverage_threshold]
         bases_above_coverage <- rowSums(bins_above_coverage)
         prop_bases_above_coverage <- bases_above_coverage / total_gene_bases
         
         data.frame(
            GeneName = rownames(depths),
            CoverageAbove = coverage_threshold,
            PropBasesAboveCoverage = prop_bases_above_coverage,
            row.names=NULL
         )
      }) 
      
      output <- do.call(rbind, output)
      
      return(output)
   }
   
   merge_tsvs(cohort_type, pattern="*.bam_metric.gene_coverage.tsv(.gz)*$", output_file_suffix = ".bam_metric.gene_coverage.tsv.gz", read_func = read_func)
}

merge_tsvs.bam_metric.exon_coverage <- function(cohort_type){ 
   merge_tsvs(cohort_type, pattern="*.bam_metric.exon_coverage.tsv(.gz)*$", output_file_suffix=".bam_metric.exon_coverage.tsv.gz") 
}

## --------------------------------
## Purple
## --------------------------------

merge_tsvs.purple.purity <- function(cohort_type){ 
   merge_tsvs(cohort_type, pattern="*.purple.purity.tsv(.gz)*$", output_file_suffix = ".purple.purity.tsv.gz") 
}

merge_tsvs.purple.qc <- function(cohort_type){
   read_func <- function(path){
      qc <- read.delim(path, header=FALSE)
      setNames(as.list(qc[,2]), qc[,1]) %>% as.data.frame()
   }
   
   merge_tsvs(cohort_type, pattern="*.purple.qc$", output_file_suffix=".purple.qc.tsv.gz", read_func=read_func)
}

merge_tsvs.purple.cnv.gene <- function(cohort_type){
   
   read_func <- function(path){
      pipe_connection <- pipe(stringr::str_glue("zcat -f {path} | cut -d '\t' -f 1-6"))
      df <- read.delim(pipe_connection, comment.char = "", colClasses = c("character", "integer", "integer", "character", "numeric", "numeric"))
      df %>% dplyr::filter(gene %in% DRIVER_GENES$CNV)
   }
   
   merge_tsvs(
      cohort_type, pattern="*.purple.cnv.gene.tsv(.gz)*$", output_file_suffix=".purple.cnv.gene.tsv.gz",
      read_func=read_func, progress_interval = 10
   )
}

## --------------------------------
## Linx
## --------------------------------

merge_tsvs.linx.driver.catalog <- function(cohort_type){ 
   merge_tsvs(cohort_type, pattern="*.linx.driver.catalog.tsv(.gz)*$", output_file_suffix=".linx.driver.catalog.tsv.gz") 
}

merge_tsvs.linx.fusion <- function(cohort_type){ 
   merge_tsvs(cohort_type, pattern="*.linx.fusion.tsv(.gz)*$", output_file_suffix=".linx.fusion.tsv.gz")
}

## --------------------------------
## Lilac
## --------------------------------

merge_tsvs.lilac.solutions <- function(cohort_type){ 
   merge_tsvs(cohort_type, pattern="*.lilac.tsv(.gz)*$", output_file_suffix=".lilac.tsv.gz") 
}

merge_tsvs.lilac.qc <- function(cohort_type){ 
   merge_tsvs(cohort_type, pattern="*.lilac.qc.tsv(.gz)*$", output_file_suffix=".lilac.qc.tsv.gz")
}

## --------------------------------
## Isofox
## --------------------------------

merge_tsvs.isofox.gene_data <- function(cohort_type){
   
   read_func <- function(path){
      sep <- if(grepl("(csv|csv.gz)$", path)){ "," } else { "\t" }
      df <- read.delim(path, sep = sep, comment.char = "", colClasses = "character")
      df %>% dplyr::select(GeneName, AdjTPM) %>% dplyr::filter(GeneName %in% DRIVER_GENES$RNA)
   }
   
   df <- merge_tsvs(
      cohort_type, pattern="*.isf.gene_data.tsv(.gz)*$", output_file_suffix=".isf.gene_data.tsv.gz",
      read_func=read_func, progress_interval=10
   )
   
   df$AdjTPM <- as.numeric(df$AdjTPM)
   
   return(df)
}

## =============================
## Load sample data
## =============================


TARGET_COHORT <- local({
   
   cohort_type <- COHORT_TYPE$TARGET
   
   tables <- list()
   tables$bam_metric.summary <- merge_tsvs.bam_metric.summary(cohort_type)
   tables$bam_metric.coverage <- merge_tsvs.bam_metric.coverage(cohort_type)
   tables$bam_metric.frag_length <- merge_tsvs.bam_metric.frag_length(cohort_type)
   tables$bam_metric.gene_coverage <- merge_tsvs.bam_metric.gene_coverage(cohort_type)
   tables$bam_metric.exon_coverage <- merge_tsvs.bam_metric.exon_coverage(cohort_type)
   tables$purple.qc <- merge_tsvs.purple.qc(cohort_type)
   tables$purple.purity <- merge_tsvs.purple.purity(cohort_type)
   tables$purple.cnv.gene <- merge_tsvs.purple.cnv.gene(cohort_type)
   tables$linx.driver.catalog <- merge_tsvs.linx.driver.catalog(cohort_type)
   tables$linx.fusion <- merge_tsvs.linx.fusion(cohort_type)
   tables$lilac.solutions <- merge_tsvs.lilac.solutions(cohort_type)
   tables$lilac.qc <- merge_tsvs.lilac.qc(cohort_type)
   tables$isofox.gene_data <- merge_tsvs.isofox.gene_data(cohort_type)
   
   return(tables)
})

REFERENCE_COHORT <- local({
   
   cohort_type <- COHORT_TYPE$REFERENCE
   
   tables <- list()
   tables$purple.purity <- merge_tsvs.purple.purity(cohort_type)
   tables$purple.cnv.gene <- merge_tsvs.purple.cnv.gene(cohort_type)
   tables$linx.driver.catalog <- merge_tsvs.linx.driver.catalog(cohort_type)
   tables$linx.fusion <- merge_tsvs.linx.fusion(cohort_type)
   tables$lilac.solutions <- merge_tsvs.lilac.solutions(cohort_type)
   tables$isofox.gene_data <- merge_tsvs.isofox.gene_data(cohort_type)
   
   return(tables)
})

## =============================
## Plots
## =============================

PLOTS_DIR <- file.path(args$output_dir, "plots")
if(!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR)

scale_values <- function(target_cohort_value, reference_cohort_value){
   values <- c(target_cohort_value, reference_cohort_value)
   names(values) <- c(COHORT_TYPE$TARGET$name, COHORT_TYPE$REFERENCE$name)
   return(values)
}

SCALE_COLOR_COHORT <- scale_color_manual(values = scale_values("#F8766D", "#00BFC4"))
SCALE_FILL_COHORT <- scale_fill_manual(values = scale_values("#F8766D", "#00BFC4"))
SCALE_SIZE_COHORT <- scale_size_manual(values = scale_values(2, 1))

SCALE_Y_LOG10 <- scale_y_log10(label = function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE))

preordered_factor <- function(x){ factor(x, unique(x)) }

merge_cohort_data <- function(target_cohort_data, reference_cohort_data){
   target_cohort_data$Cohort <- COHORT_TYPE$TARGET$name
   reference_cohort_data$Cohort <- COHORT_TYPE$REFERENCE$name
   
   df <- dplyr::bind_rows(target_cohort_data, reference_cohort_data)
   df$Cohort <- preordered_factor(df$Cohort)
   
   return(df)
}

form_plot_path <- function(output_file_suffix){
   file.path(PLOTS_DIR, paste0(COHORT_TYPE$TARGET$name , output_file_suffix))
}

plot_to_pdf <- function(plots, output_path, width, height){
   
   LOGGER$info("Writing PDF: %s", output_path)
   
   pdf(output_path, width = width, height = height)
   
   if(is.list(plots)){
      for(i in seq_along(plots)){
         LOGGER$debug("  Writing page: %s", i)
         page <- plots[[i]]
         plot(page)
      }
   } else {
      plot(plots)
   }
   
   dev.off()
}

## --------------------------------
## Summary
## --------------------------------

percentile_rank <- function(x){
   #x = c(2,1,0,0,3,4,5)
   (rank(x, ties="first")-1) / (length(x)-1)
}

plot_ecdf <- function(plot_data, y, color = NULL){
   
   if(is.null(color)){
      plot_data <- plot_data %>% dplyr::mutate(Rank = percentile_rank(.data[[y]]))
   } else {
      plot_data <- plot_data %>% dplyr::group_by(.data[[color]]) %>% dplyr::mutate(Rank = percentile_rank(.data[[y]]))
   }
   
   plot <- 
      ggplot(plot_data, aes(
         x = Rank, y = .data[[y]],
         color = if(!is.null(color)){ .data[[color]] } else { NULL }
      )) +
      geom_line()
   
   if(!is.null(color)){
      plot <- plot + labs(color = color)
   }
   
   return(plot)
}

SUMMARY_PLOTS <- local({
   
   plots <- list()

   ## Row 1
   plots$mean_coverage <- TARGET_COHORT$bam_metric.summary %>%
      plot_ecdf(y = "MeanCoverage") +
      coord_cartesian(ylim = c(0, NA))

   plots$total_reads <- TARGET_COHORT$bam_metric.summary %>%
      plot_ecdf(y = "TotalReads") +
      coord_cartesian(ylim = c(0, NA))

   plots$on_target_rate <- TARGET_COHORT$bam_metric.summary %>%
      dplyr::mutate(PropOnTargetReads = 1 - OffTargetReads / (TotalReads+OffTargetReads) ) %>%
      plot_ecdf(y = "PropOnTargetReads") +
      coord_cartesian(ylim = c(0, 1))

   plots$duplicate_rate <- TARGET_COHORT$bam_metric.summary %>%
      dplyr::mutate(PropDuplicateReads = DuplicateReads / TotalReads) %>%
      plot_ecdf(y = "PropDuplicateReads") +
      coord_cartesian(ylim = c(0, 1))

   plots$fragment_length <- TARGET_COHORT$bam_metric.frag_length %>%
      plot_ecdf(y = "MedianFragmentLength") +
      coord_cartesian(ylim = c(0, NA))
   
   plots$coverage_above <- TARGET_COHORT$bam_metric.coverage %>%
      dplyr::mutate(CoverageAbove = factor(CoverageAbove, sort(unique(CoverageAbove)))) %>%
      plot_ecdf(y = "PropBasesAboveCoverage", color = "CoverageAbove") +
      coord_cartesian(ylim = c(0, 1))
   
   purple_purity <- merge_cohort_data(TARGET_COHORT$purple.purity, REFERENCE_COHORT$purple.purity)
   
   ## Row 2
   plots$purity <- purple_purity %>% 
      dplyr::filter(status != "NO_TUMOR") %>%
      plot_ecdf(y = "purity", color = "Cohort") + 
      coord_cartesian(ylim = c(0, 1))
   
   plots$ploidy <- purple_purity %>% 
      dplyr::filter(status != "NO_TUMOR") %>%
      plot_ecdf(y = "ploidy", color = "Cohort") + 
      coord_cartesian(ylim = c(0, NA))
   
   plots$tmb <- purple_purity %>% 
      dplyr::filter(status != "NO_TUMOR") %>%
      plot_ecdf(y = "tmbPerMb", color = "Cohort") + 
      SCALE_Y_LOG10
   
   plots$ms_indels_per_mb <- purple_purity %>% 
      dplyr::filter(status != "NO_TUMOR") %>%
      plot_ecdf(y = "msIndelsPerMb", color = "Cohort") + 
      SCALE_Y_LOG10
   
   ## Row 3
   plots$qc_status <- TARGET_COHORT$purple.qc %>% 
      dplyr::count(QCStatus, name = "Count") %>%
      ggplot(aes(y = QCStatus, x = Count)) +
      geom_bar(stat = "identity", fill = "grey") +
      geom_text(aes(x = 0, label = paste0(Count, ": ", QCStatus)), hjust = 0, size = 2) +
      labs(y = "PurpleQcStatus") +
      theme(
         panel.grid.major.y = element_blank(),
         panel.grid.minor.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
      )
   
   plots$no_tumor <- purple_purity %>% 
      dplyr::group_by(Cohort) %>% 
      dplyr::summarise(PropNoTumor = sum(status == "NO_TUMOR") / length(SampleId)) %>%
      ggplot(aes(x = Cohort, y = PropNoTumor, fill = Cohort)) + 
      geom_bar(stat="identity") + 
      geom_text(aes(y = 0, label = round(PropNoTumor, 3)), vjust = 0, size = 2.7) +
      coord_cartesian(ylim = c(0, NA)) +
      theme(legend.position = "none")
   
   plots$qc_status_lilac <- TARGET_COHORT$lilac.qc %>% 
      dplyr::rename(QCStatus = "Status") %>%
      dplyr::count(QCStatus, name = "Count") %>%
      ggplot(aes(y = QCStatus, x = Count)) +
      geom_bar(stat = "identity", fill = "grey") +
      geom_text(aes(x = 0, label = paste0(Count, ": ", QCStatus)), hjust = 0, size = 2) +
      labs(y = "LilacQcStatus") +
      theme(
         panel.grid.major.y = element_blank(),
         panel.grid.minor.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
      )

   patchwork::wrap_plots(plots, ncol = 5, guides = "collect")
})

plot_to_pdf(SUMMARY_PLOTS, form_plot_path(".qc_summary.pdf"), width = 15, height = 7)

## --------------------------------
## Exon coverage
## --------------------------------

exon_key <- function(df){ paste(df$GeneName, df$ExonRank) }

split_genes_over_pages <- function(genes, rows = 8, cols = 12) {
   #genes <- unique(TARGET_COHORT$bam_metric.gene_coverage$GeneName)
   
   max_genes_per_page <- rows * cols
   
   gene_page_numbers <- ceiling(seq_along(genes) / max_genes_per_page)
   genes_per_page <- split(genes, gene_page_numbers)
   
   config_per_page <- lapply(genes_per_page, function(page_genes){
      rows_to_pad <- 0
      
      n_page_genes <- length(page_genes)
      if(n_page_genes < max_genes_per_page){
         n_missing <- max_genes_per_page - n_page_genes
         rows_to_pad <- floor(n_missing / cols)
      }
      
      return(list(
         genes = page_genes,
         cols = cols,
         rows = rows - rows_to_pad,
         rows_to_pad = rows_to_pad
      ))
   })
   
   return(config_per_page)
}


EXON_INFO <- local({
   
   ## Extract exon gene and coords. Could also get this info from ensembl data cache
   exons <- TARGET_COHORT$bam_metric.exon_coverage %>%
      dplyr::select(GeneName, Chromosome, PosStart, PosEnd, ExonRank) %>%
      unique()
   
   exons <- exons[naturalsort::naturalorder(exons$Chromosome),]
   
   ## Convert to genomic ranges
   exons_gr <- GenomicRanges::GRanges(
      seqnames = exons$Chromosome,
      ranges = IRanges::IRanges(start = exons$PosStart, end = exons$PosEnd),
      gene = exons$GeneName
   )
   
   panel_gr <- GenomicRanges::GRanges(
      seqnames = PANEL_BED$Chromosome,
      ranges = IRanges::IRanges(start = PANEL_BED$PosStart, end = PANEL_BED$PosEnd),
      description = PANEL_BED$description
   )
   
   ## Match exons in samples to panel regions and get coverage
   hits <- GenomicRanges::findOverlaps(query = exons_gr, subject = panel_gr)
   exon_hits <- S4Vectors::queryHits(hits)
   panel_hits <- S4Vectors::subjectHits(hits)
   
   hits_coverage <- cbind(
      exons[exon_hits,],
      PANEL_BED[panel_hits,] %>% dplyr::select(-Chromosome) %>% dplyr::rename_with(~ paste0("Panel", .))
   )
   
   hits_coverage <- hits_coverage %>% dplyr::mutate(
      MissingStart = pmax(PanelPosStart - PosStart, 0),
      MissingEnd = pmax(PosEnd - PanelPosEnd, 0),
      TotalBases = PosEnd - PosStart,
      CoveredBases = TotalBases - MissingStart - MissingEnd,
      CoveredProp = CoveredBases / TotalBases
   )
   
   ## Add 0 coverage for non hits
   non_hits_coverage <- exons[!(exon_key(exons) %in% exon_key(hits_coverage)),]
   non_hits_coverage$CoveredProp <- 0
   exon_coverage <- dplyr::bind_rows(hits_coverage, non_hits_coverage)
   
   ## Output
   exon_coverage %>% dplyr::arrange(GeneName, ExonRank)
})

EXON_COVERAGE_PLOTS <- local({

   plot_data <- TARGET_COHORT$bam_metric.exon_coverage %>%
      dplyr::group_by(GeneName, ExonRank) %>%
      dplyr::summarise(MeanPropExonAbove250X = mean(PercAboveDepth_250)) %>%
      dplyr::ungroup()
   
   plot_data$ExonCoveredProp <- EXON_INFO$CoveredProp[match(exon_key(plot_data), exon_key(EXON_INFO))]
   plot_data$ExonCoveredPerc <- 100 * round(plot_data$ExonCoveredProp, 2)
   plot_data$ExonCoveredPerc[plot_data$ExonCoveredProp == 0] <- ""
   
   GENES_PER_PAGE <- 100
   page_configs <- plot_data$GeneName %>% unique() %>% sort() %>% split_genes_over_pages(rows=GENES_PER_PAGE, cols=1)

   pages <- lapply(page_configs, function(page_config){
      #page_config <- page_configs[[length(page_configs)]]
      
      plot <- plot_data %>%
         
         dplyr::filter(GeneName %in% page_config$genes) %>%

         ggplot(aes(x = ExonRank, y = GeneName)) +
         coord_cartesian(expand = FALSE) +

         geom_tile(aes(fill = MeanPropExonAbove250X), color = "black", linewidth = 0.1) +
         geom_text(aes(label = ExonCoveredPerc), size = 2) +

         scale_y_discrete(limits = rev) +
         scale_x_continuous(
            breaks = seq(10, 1000, by = 10),
            limits = c(0.5, max(plot_data$ExonRank+0.5))
         ) +
         scale_fill_distiller(palette = "Spectral", limit = c(0, NA), na.value = "white") +

         labs(title = "Label overlay: % of exon_coverage.tsv regions covered by panel bed file") +

         theme(
            plot.title = element_text(size = 11),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            legend.position = "top",
            legend.justification = "left"
         )

      patchwork::wrap_plots(plot, patchwork::plot_spacer(), heights = c(page_config$rows, page_config$rows_to_pad))
   })

   attr(pages, "width") <- max(plot_data$ExonRank)
   attr(pages, "height") <- GENES_PER_PAGE

   return(pages)
})

plot_to_pdf(
   EXON_COVERAGE_PLOTS,
   form_plot_path(".exon_coverage.pdf"),
   width = attr(EXON_COVERAGE_PLOTS, "width") * 0.2,
   height = attr(EXON_COVERAGE_PLOTS, "height") * 0.2
)


## --------------------------------
## Gene coverage
## --------------------------------

GENE_COVERAGE_PLOTS <- local({
   
   page_configs <- TARGET_COHORT$bam_metric.gene_coverage$GeneName %>% unique() %>% split_genes_over_pages(rows = 8, cols = 12)
   
   pages <- lapply(page_configs, function(page_config){
      #page_config=page_configs[[length(page_configs)]]
      plot_data <- TARGET_COHORT$bam_metric.gene_coverage %>% 
         dplyr::filter(GeneName %in% page_config$genes) %>%
         dplyr::group_by(GeneName, CoverageAbove) %>%
         dplyr::mutate(Rank = percentile_rank(PropBasesAboveCoverage)) %>%
         dplyr::ungroup() %>%
         dplyr::mutate(CoverageAbove = factor(CoverageAbove, sort(unique(CoverageAbove))))
      
      plot <- ggplot(plot_data, aes(x = Rank, y = PropBasesAboveCoverage, color = CoverageAbove)) +
         geom_line() +
         facet_wrap(GeneName ~ ., ncol = page_config$cols) +
         coord_cartesian(ylim = c(0, 1)) +
         theme(
            panel.spacing = unit(1, "pt"),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
         )
      
      patchwork::wrap_plots(plot, patchwork::plot_spacer(), heights = c(page_config$rows, page_config$rows_to_pad))
   })
   
   return(pages)
})

plot_to_pdf(GENE_COVERAGE_PLOTS, form_plot_path(".gene_coverage.pdf"), height = 8, width = 12)


## --------------------------------
## Driver
## --------------------------------

DRIVER_FREQ_PLOT <- local({
   
   MIN_DRIVER_LIKELIHOOD <- 0.2
   SELECTED_DRIVER_TYPES <- c("MUTATION", "AMP", "DEL")
   
   calc_driver_freq <- function(driver_catalog, n_samples){
      #driver_catalog=TARGET_COHORT$linx.driver.catalog
      #n_samples=nrow(TARGET_COHORT$purple.purity)
      
      driver_catalog <- driver_catalog %>% dplyr::filter(
         driverLikelihood > MIN_DRIVER_LIKELIHOOD & 
         driver %in% SELECTED_DRIVER_TYPES &
         isCanonical == "true"
      )
      
      if("reportedStatus" %in% colnames(driver_catalog)){
         driver_catalog <- driver_catalog %>% dplyr::filter(reportedStatus == "REPORTED")
      }
      
      driver_freq <- driver_catalog %>% 
         dplyr::group_by(gene, driver) %>% 
         dplyr::summarise(Frequency = sum(driverLikelihood) / n_samples, .groups = "drop_last")
      
      return(driver_freq)
   }
   
   driver_freqs <- list()
   driver_freqs$target <- calc_driver_freq(TARGET_COHORT$linx.driver.catalog, nrow(TARGET_COHORT$purple.purity))
   driver_freqs$reference <- calc_driver_freq(REFERENCE_COHORT$linx.driver.catalog, nrow(REFERENCE_COHORT$purple.purity))
   
   ## Use merge to:
   ## - Select only the driver events present in the target cohort
   ## - Fill in driver events missing in reference cohort
   driver_freq_wide <- merge(driver_freqs$target , driver_freqs$reference, by = c("gene", "driver"), all.x = TRUE)
   colnames(driver_freq_wide) <- c(c("gene", "driver", COHORT_TYPE$TARGET$name, COHORT_TYPE$REFERENCE$name))
   driver_freq_wide[[COHORT_TYPE$REFERENCE$name]][ is.na(driver_freq_wide[[COHORT_TYPE$REFERENCE$name]]) ] <- 0
   
   ## Plot
   driver_freq_long <- driver_freq_wide %>% reshape2::melt(
      id.vars = c("gene", "driver"), 
      variable.name = "Cohort", 
      value.name = "Frequency"
   )
   
   plots <- lapply(SELECTED_DRIVER_TYPES, function(driver_type){
      #driver_type="MUTATION"
      plot_data_wide <- driver_freq_wide %>% dplyr::filter(driver == driver_type)
      plot_data_long <- driver_freq_long %>% dplyr::filter(driver == driver_type)
      
      gene_order <- plot_data_wide[order(plot_data_wide[[COHORT_TYPE$TARGET$name]], decreasing = TRUE),"gene"]
      plot_data_long$gene <- factor(plot_data_long$gene, gene_order)
      
      ggplot(plot_data_long, aes(x=gene, y=Frequency)) + 
         geom_point(aes(color=Cohort, size=Cohort)) + 
         SCALE_COLOR_COHORT +
         SCALE_SIZE_COHORT +
         labs(title = driver_type) +
         theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
         )
   })
   
   patchwork::wrap_plots(plots, ncol = 1, guides = "collect", axes = "collect")
})

plot_to_pdf(DRIVER_FREQ_PLOT, form_plot_path(".frequency.driver.pdf"), width = 25, height = 7)

## --------------------------------
## Fusions
## --------------------------------

if(!is.null(TARGET_COHORT$linx.fusion)){

   REPORTABLE_FUSIONS_PLOT <- local({
      
      if(is.null(TARGET_COHORT$linx.fusion))
         return(NULL)
      
      calc_fusion_freq <- function(fusions, n_samples){
         fusions %>% 
            dplyr::filter(reported=="true") %>% 
            dplyr::count(name, name = "Count") %>% 
            dplyr::mutate(Frequency = Count / n_samples)
      }
      
      fusion_freqs <- list()
      fusion_freqs$target <- calc_fusion_freq(TARGET_COHORT$linx.fusion, nrow(TARGET_COHORT$purple.purity))
      fusion_freqs$reference <- calc_fusion_freq(REFERENCE_COHORT$linx.fusion, nrow(REFERENCE_COHORT$purple.purity))
      fusion_freqs$reference <- fusion_freqs$reference %>% dplyr::filter(name %in% fusion_freqs$target$name)
      
      fusion_freq <- merge_cohort_data(fusion_freqs$target, fusion_freqs$reference)
      
      fusions_order <- fusion_freqs$target %>% dplyr::arrange(Frequency) %>% pull(name)
      fusion_freq <- fusion_freq %>% dplyr::mutate(
         name = factor(name, fusions_order),
         Cohort = preordered_factor(Cohort)
      )
      
      plot <- ggplot(fusion_freq, aes(x = name, y = Frequency)) +
         geom_point(aes(color = Cohort, size = Cohort)) +
         SCALE_SIZE_COHORT +
         coord_cartesian(ylim = c(0, NA)) +
         theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
         )
      
      plot$width <- length(fusions_order)
      
      return(plot)
   })
      

   plot_to_pdf(
      REPORTABLE_FUSIONS_PLOT, 
      form_plot_path(".frequency.fusion.pdf"), 
      width = REPORTABLE_FUSIONS_PLOT$width*0.5 + 1, 
      height = 5
   )
   
}

## --------------------------------
## Gene copy number
## --------------------------------

plot_gene_cnv <- function(cnv_type){
   
   if(cnv_type == "AMP"){
      selected_genes <- DRIVER_GENES$AMP
      y_variable <- "relativeMaxCopyNumber"
   } else if(cnv_type == "DEL"){
      selected_genes <- DRIVER_GENES$DEL
      y_variable <- "minCopyNumber"
   } else {
      stop("Invalid cnv_type: ", cnv_type)
   }
   
   select_top_n_samples <- function(df, top_n){
      set.seed(0)
      selected_samples <- df$SampleId %>% unique() %>% sample(top_n)
      df %>% dplyr::filter(SampleId %in% selected_samples)
   }
   
   cnv_data <- list()
   
   cnv_data$target <- TARGET_COHORT$purple.cnv.gene
   cnv_data$reference <- REFERENCE_COHORT$purple.cnv.gene %>% select_top_n_samples(1000)
   
   cnv_data$target$ploidy <- TARGET_COHORT$purple.purity$ploidy[ match(cnv_data$target$SampleId, TARGET_COHORT$purple.purity$SampleId) ]
   cnv_data$reference$ploidy <- TARGET_COHORT$purple.purity$ploidy[ match(cnv_data$reference$SampleId, REFERENCE_COHORT$purple.purity$SampleId) ]
   
   cnv_data <- merge_cohort_data(cnv_data$target, cnv_data$reference)
   
   cnv_data <- cnv_data %>% 
      dplyr::filter(gene %in% selected_genes) %>%
      dplyr::mutate(
         relativeMaxCopyNumber = maxCopyNumber / ploidy,
         relativeMaxCopyNumber = pmax(0, relativeMaxCopyNumber)
      ) %>%
      dplyr::mutate(
         Cohort = preordered_factor(Cohort),
         gene = preordered_factor(gene)
      ) %>%
      dplyr::mutate(
         chromosome = ifelse(grepl("^chr", chromosome), chromosome, paste0("chr", chromosome)),
         chromosome = preordered_factor(chromosome)
      )
   
   chromosomes <- levels(cnv_data$chromosome)
   pages <- lapply(chromosomes, function(chrom){
      #chrom <- "chr20"
      plot <- cnv_data %>% 
         dplyr::filter(chromosome == chrom) %>%
         ggplot(aes(x = gene, y = .data[[y_variable]], fill = Cohort)) + 
         geom_jitter(position = position_jitterdodge(jitter.width = 0.2, seed = 0), size = 0.2, alpha = 0.3) +
         geom_boxplot(outlier.color = NA, alpha = 0.5, linewidth = 0.2, show.legend = TRUE) +
         scale_fill_discrete(drop = FALSE) +
         labs(title = chrom) +
         theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
         )
      
      if(cnv_type == "AMP"){
         plot <- plot +
            SCALE_Y_LOG10 +
            geom_hline(yintercept = 3, linetype="11", color = "red") +
            geom_hline(yintercept = 3*0.8, linetype="11", color = "orange") +
            coord_cartesian(ylim = c(0.01, 100))
      }
      
      if(cnv_type == "DEL"){
         plot <- plot + 
            geom_hline(yintercept = 0.5, linetype="11", color = "red") +
            coord_cartesian(ylim = c(0, 6))
      }
      
      return(plot)
   })
   
   return(pages)
   
}

GENE_AMP_PLOTS <- plot_gene_cnv("AMP")
GENE_DEL_PLOTS <- plot_gene_cnv("DEL")

plot_to_pdf(GENE_AMP_PLOTS, form_plot_path(".gene_cnv.amp.pdf"), width = 10, height = 3)
plot_to_pdf(GENE_DEL_PLOTS, form_plot_path(".gene_cnv.del.pdf"), width = 10, height = 3)

## --------------------------------
## Lilac
## --------------------------------

if(!is.null(TARGET_COHORT$lilac.solutions)){
   
   HLA_ALLELE_PLOT <- local({
      
      calc_allele_freq <- function(df){ df %>% dplyr::count(Allele, name = "Count") %>% dplyr::mutate(CohortFrequency = Count/sum(Count)) }
      
      allele_freqs <- list()
      
      allele_freqs$target <- calc_allele_freq(TARGET_COHORT$lilac.solutions)
      allele_freqs$reference <- calc_allele_freq(REFERENCE_COHORT$lilac.solutions)
      allele_freqs$reference <- allele_freqs$reference %>% dplyr::filter(Allele %in% allele_freqs$target$Allele)
      
      allele_freq <- merge_cohort_data(allele_freqs$target, allele_freqs$reference)
      
      allele_order <- allele_freqs$target %>% dplyr::arrange(CohortFrequency) %>% pull(Allele)
      allele_freq <- allele_freq %>% dplyr::mutate(Allele = factor(Allele, allele_order))
      
      qc_status_counts <- TARGET_COHORT$lilac.qc %>% dplyr::count(Status)
      qc_status_string <- paste0(qc_status_counts$Status,"=",qc_status_counts$n) %>% paste(collapse = '\n')
      
      plot <- ggplot(allele_freq, aes(x = Allele, y = CohortFrequency)) +
         geom_point(aes(color = Cohort, size = Cohort)) +
         SCALE_COLOR_COHORT +
         SCALE_SIZE_COHORT +
         coord_cartesian(ylim = c(0, NA)) +
         labs(title = qc_status_string) +
         theme(
            plot.title = element_text(size = 8),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
         )
      
      plot$width <- allele_freq$Allele %>% unique() %>% length()
      
      return(plot)
   })
   
   plot_to_pdf(HLA_ALLELE_PLOT, form_plot_path(".frequency.hla_allele.pdf"), width = HLA_ALLELE_PLOT$width*0.15, height = 6)
}

## --------------------------------
## Plots - RNA (TPM)
## --------------------------------

if(!is.null(TARGET_COHORT$isofox.gene_data)){
   
   RNA_ADJ_TPM_PLOTS <- local({
      
      gene_data <- merge_cohort_data(
         TARGET_COHORT$isofox.gene_data,
         REFERENCE_COHORT$isofox.gene_data %>% dplyr::filter(GeneName %in% TARGET_COHORT$isofox.gene_data$GeneName)
      )
      
      page_configs <- gene_data$GeneName %>% unique() %>% split_genes_over_pages(rows = 8, cols = 12)
      
      pages <- lapply(page_configs, function(page_config){
         #page_config=page_configs[[1]]
         plot_data <- gene_data %>%
            dplyr::filter(GeneName %in% page_config$genes) %>%
            dplyr::group_by(GeneName, Cohort) %>%
            dplyr::mutate(Rank = percentile_rank(AdjTPM)) %>%
            dplyr::ungroup()
         
         plot <- ggplot(plot_data, aes(x = Rank, y = log10(AdjTPM + 0.01), color = Cohort)) +
            geom_line() +
            facet_wrap(GeneName ~ ., ncol = page_config$cols) +
            theme(
               panel.spacing = unit(1, "pt"),
               axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            )
         
         patchwork::wrap_plots(plot, patchwork::plot_spacer(), heights = c(page_config$rows, page_config$rows_to_pad))
      })
      
      return(pages)
   })
      
   plot_to_pdf(RNA_ADJ_TPM_PLOTS, form_plot_path(".rna.adj_tpm.pdf"), width = 12, height = 8)
}
