options(warn = 1)

if(!interactive()){
   ## Prevent empty file Rplots.pdf from being written
   ## See: https://stackoverflow.com/questions/6535927/how-do-i-prevent-rplots-pdf-from-being-generated
   pdf(NULL)
}

suppressPackageStartupMessages(library(dplyr))

library(ggplot2)
theme_set(
   theme_bw() +
   theme(panel.grid = element_blank())
)

library(patchwork)
library(scales)

## =============================
## Config
## =============================

args <- commandArgs(trailingOnly = TRUE)

parseOptionalArg <- function(arg){
   if(arg == "NA") NA  else arg
}

TUMOR_ID <- args[1]
NORMAL_ID <- parseOptionalArg(args[2])
VIS_DATA_FILE <- args[3]
COHORT_PERCENTILES_FILE <- parseOptionalArg(args[4])
PLOT_PATH <- args[5]
SHOW_PLOT_WARNINGS <- as.logical(args[6])
SINGLE_PATIENT_MODE <- as.logical(args[7])
GLOBAL_LOG_LEVEL <- args[8]

if(FALSE){
   TUMOR_ID <- "TUMOR"
   NORMAL_ID <- "TUMOR-ref"
   VIS_DATA_FILE <- "TUMOR.qsee.vis.data.tsv.gz"
   COHORT_PERCENTILES_FILE <- "qsee.cohort.percentiles.tsv.gz"
   PLOT_PATH <- "TUMOR.qsee.vis.report.pdf"
   
   SHOW_PLOT_WARNINGS <- TRUE
   SINGLE_PATIENT_MODE <- TRUE
   GLOBAL_LOG_LEVEL <- "DEBUG"
}

## =============================
## Logging
## =============================

GLOBAL_LOG_LEVEL <- toupper(GLOBAL_LOG_LEVEL)

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

   log_prefix <- if(SINGLE_PATIENT_MODE) "" else sprintf("[%s] ", TUMOR_ID)

   log_message <- sprintf("%s [R] [%-5s] %s%s",
      current_time, log_level$name, log_prefix, sprintf(fmt, ...))
   
   if(log_level$severity >= LOG_LEVEL[[LOG_LEVEL$ERROR$name]]$severity)
      stop(log_message)
   
   if(log_level$severity >= LOG_LEVEL[[GLOBAL_LOG_LEVEL]]$severity)
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

print_arg <- function(arg){
   LOGGER$debug("  %s: %s", tolower(deparse(substitute(arg))), arg)
}

LOGGER$debug("Running script with args:")
print_arg(TUMOR_ID)
print_arg(NORMAL_ID)
print_arg(VIS_DATA_FILE)
print_arg(COHORT_PERCENTILES_FILE)
print_arg(PLOT_PATH)
print_arg(SHOW_PLOT_WARNINGS)
print_arg(SINGLE_PATIENT_MODE)
print_arg(GLOBAL_LOG_LEVEL)

## =============================
## Constants
## =============================

SAMPLE_TYPE <- list(
   TUMOR = list(name = "TUMOR", display_name = "Tumor", color = "#D1392C"),
   NORMAL = list(name = "NORMAL", display_name = "Normal", color = "#4A7DB4")
)

FEATURE_TYPE <- list(
   SUMMARY_TABLE              = "SUMMARY_TABLE",
   COVERAGE_DISTRIBUTION      = "COVERAGE_DISTRIBUTION",
   FRAG_LENGTH_DISTRIBUTION   = "FRAG_LENGTH_DISTRIBUTION",
   GC_BIAS                    = "GC_BIAS",
   DUPLICATE_FREQ             = "DUPLICATE_FREQ",
   DISCORDANT_FRAG_FREQ       = "DISCORDANT_FRAG_FREQ",
   MISSED_VARIANT_LIKELIHOOD  = "MISSED_VARIANT_LIKELIHOOD",
   BQR_BY_SNV96_CONTEXT       = "BQR_PER_SNV96_CONTEXT",
   BQR_BY_ORIG_QUAL           = "BQR_PER_ORIG_QUAL",
   MS_INDEL_ERROR_RATES       = "MS_INDEL_ERROR_RATES",
   MS_INDEL_ERROR_BIAS        = "MS_INDEL_ERROR_BIAS"
)

PERCENTILE_PREFIX <- "Pct_"

NAMED_PERCENTILES <- list(
   MIN   = list(name = "PctMin"  , colname = paste0(PERCENTILE_PREFIX,  5)),
   LOWER = list(name = "PctLower", colname = paste0(PERCENTILE_PREFIX, 25)),
   MID   = list(name = "PctMid"  , colname = paste0(PERCENTILE_PREFIX, 50)),
   UPPER = list(name = "PctUpper", colname = paste0(PERCENTILE_PREFIX, 75)),
   MAX   = list(name = "PctMax"  , colname = paste0(PERCENTILE_PREFIX, 95))
)

QC_STATUS <- list(
   FAIL = list(name = "FAIL", color = "#F3C27B"),
   WARN = list(name = "WARN", color = "#FCFFC6"),
   PASS = list(name = "PASS", color = "#FFFFFF00")
)

## =============================
## Load data
## =============================

load_cohort_percentiles <- function(){
   
   if(is.na(COHORT_PERCENTILES_FILE)){
      LOGGER$info("No cohort percentiles provided - skipping plotting cohort data")
      return(NULL)
   }
   
   LOGGER$info("Loading cohort percentiles from: %s", COHORT_PERCENTILES_FILE)
   
   cohort_percentiles <- read.delim(COHORT_PERCENTILES_FILE, na.strings = c("NA", "null"))
   
   cohort_named_percentiles <- cohort_percentiles[ sapply(NAMED_PERCENTILES, `[[`, "colname") ]
   colnames(cohort_named_percentiles) <- sapply(NAMED_PERCENTILES, `[[`, "name")
   
   cohort_named_percentiles <- cbind(
      cohort_percentiles %>% dplyr::select(!dplyr::starts_with(PERCENTILE_PREFIX)),
      cohort_named_percentiles
   )
   
   return(cohort_named_percentiles)
}

load_sample_features <- function(){
   LOGGER$info("Loading sample features from: %s", VIS_DATA_FILE)
   
   sample_features <- read.delim(VIS_DATA_FILE, na.strings = c("NA", "null"))
   sample_features <- sample_features %>% dplyr::filter(SampleId %in% c(TUMOR_ID, NORMAL_ID))
   
   return(sample_features)
}

COHORT_DATA <- load_cohort_percentiles()
SAMPLE_DATA <- load_sample_features()

PLOTS <- list()

## =============================
## String functions
## =============================

KEY_VALUE_SEP <- ":"
GROUP_SEP <- ";"

df_to_strings <- function(df, key_value_sep = KEY_VALUE_SEP, group_sep = GROUP_SEP){
   
   key_value_strings <- lapply(colnames(df), function(colname){
      paste0(colname, key_value_sep, df[[colname]])
   })
   
   strings <- do.call(paste, c(key_value_strings, sep=group_sep))
   
   return(strings)
}

strings_to_df <- function(strings, key_value_sep = KEY_VALUE_SEP, group_sep = GROUP_SEP){
   
   strings_split <- strsplit(as.character(strings), group_sep, fixed=TRUE)
   
   key_value_pairs <- lapply(strings_split, function(row_strings){
      row_strings_split <- strsplit(row_strings, key_value_sep)
      keys <- sapply(row_strings_split, `[[`, 1)
      values <- sapply(row_strings_split, `[[`, 2)
      structure(values, names = keys)
   })
   
   df <- dplyr::bind_rows(key_value_pairs) %>% as.data.frame(check.names=FALSE)
   df[is.na(df)] <- ""
   
   return(df)
}

has_multiple_fields <- function(strings, key_value_sep = KEY_VALUE_SEP, group_sep = GROUP_SEP){
   
   multifield_string_regex <- paste0(".+", key_value_sep, ".+", group_sep, "*")
   is_multifield_string <- grepl(multifield_string_regex, strings)
   
   if(unique(is_multifield_string) > 1){
      LOGGER$error("Found a mix of single and multi field strings: %s", paste(strings, collapse = "\n"))
   }
   
   return(all(is_multifield_string))
}

preordered_factor <- function(x){ factor(x, unique(x)) }

preordered_factors <- function(df){ 
   as.data.frame(lapply(df, function(x){
      if(is.character(x)){
         x <- factor(x, unique(x))
      }
      return(x)
   }))
}

reverse_levels <- function(fct){
   factor(fct, rev(levels(fct)))
}

## =============================
## Plot helper functions
## =============================

get_plot_data <- function(feature_type){
   
   LOGGER$info("Plotting featureType(%s)", feature_type)
   
   ## Select and merge sample/cohort data
   sample_data <- SAMPLE_DATA %>% dplyr::filter(FeatureType == feature_type) %>% select(-SourceTool)
   
   if(nrow(sample_data) == 0){
      return(NULL)
   }
   
   if(!is.null(COHORT_DATA)){
      cohort_data <- COHORT_DATA %>% dplyr::filter(FeatureType == feature_type) %>% select(-SourceTool)
      
      plot_data <- merge(
         sample_data, cohort_data, 
         by=c("SampleType", "FeatureType", "FeatureName"), 
         all.x = TRUE, ## Only select the features in the cohort data that are present in the sample
         sort = FALSE
      )
   } else {
      plot_data <- sample_data
      
      for(named_percentile in NAMED_PERCENTILES){
         plot_data[[named_percentile$name]] <- NA
      }
   }

   ## Split multi-field strings
   for(colname in colnames(plot_data)){
      column <- plot_data[[colname]]
      
      if(!has_multiple_fields(column))
         next
      
      column_as_df <- strings_to_df(column)
      
      plot_data[[colname]] <- NULL
      plot_data <- cbind(plot_data, column_as_df)
   }

   plot_data <- preordered_factors(plot_data)
   
   return(plot_data)
}

plot_missing_data <- function(plot_labels = labs()){
   
   LOGGER$warn("  Missing sample data - creating empty plot")
   
   ggplot() +
      annotate("text", x = 0, y = 0, label = "Missing sample data") +
      plot_labels +
      theme(
         axis.text = element_blank(),
         axis.ticks = element_blank()
      )
}

render_now <- function() structure(list(), class = "render_now")

ggplot_add.render_now <- function(object, plot, object_name) {
   #' Force warnings to be shown immediately.
   #' Usage: plot + render_now()
   if(SHOW_PLOT_WARNINGS){
      ggplot2::ggplotGrob(plot)
   }

   plot
}

## =============================
## Line / PDF
## =============================

plot_distribution <- function(plot_data, x, invert_normal = FALSE, mark_sample_peak = FALSE, hlines = NULL){
   
   if(FALSE){
      plot_data = get_plot_data(FEATURE_TYPE$COVERAGE_DISTRIBUTION)
      x = "ReadDepth"
   }
   
   if(is.null(plot_data)){
      LOGGER$error("Line/distribution plot failed - empty data frame")
   }
   
   plot_data[[x]] <- as.numeric(as.character(plot_data[[x]]))
   
   gg_scale_y_continuous <- geom_blank()
   sample_type_count <- plot_data$SampleType %>% levels() %>% length()
   if(invert_normal && sample_type_count > 1){
      
      signs <- ifelse(plot_data$SampleType == SAMPLE_TYPE$NORMAL$name, -1, 1)
      
      plot_data <- plot_data %>% dplyr::mutate(
         PctMin   = signs * PctMin,
         PctLower = signs * PctLower,
         PctMid   = signs * PctMid,
         PctUpper = signs * PctUpper,
         PctMax   = signs * PctMax,
         FeatureValue = signs * FeatureValue
      )
      
      gg_scale_y_continuous <- scale_y_continuous(labels = function(x) abs(x))
   }
   
   gg_geom_line <- geom_line(aes(color = SampleType))
   gg_geom_segment <- geom_blank()
   if(mark_sample_peak){
      
      peak_data <- plot_data %>% 
         dplyr::group_by(SampleType) %>%
         
         ## Remove the bottom 5% since the coverage plot has a spike there
         dplyr::slice(floor(dplyr::n() * 0.05):ceiling(dplyr::n() * 0.95)) %>% 
         
         dplyr::reframe(
            XPos = .data[[x]][which.max(abs(FeatureValue))],
            Height = FeatureValue[which.max(abs(FeatureValue))]
         )
      
      gg_geom_segment <- geom_segment(
         data = peak_data, 
         mapping = aes(x = XPos, xend = XPos, y = 0, yend = Height, color = SampleType),
         show.legend = FALSE
      )
   }
   
   gg_geom_ribbon <- geom_ribbon(aes(ymin = PctMin, ymax = PctMax, fill = SampleType), alpha = 0.1)
   
   sample_type_colors <- sapply(SAMPLE_TYPE, `[[`, "color")
   gg_scale_fill_manual <- scale_fill_manual(values = sample_type_colors)
   gg_scale_color_manual <- scale_color_manual(values = sample_type_colors)
   
   if(is.na(COHORT_PERCENTILES_FILE)){
      gg_geom_ribbon <- geom_blank()
      gg_scale_fill_manual <- geom_blank()
   }
   
   ggplot(plot_data, aes(x = .data[[x]], y = FeatureValue, group = SampleType)) +
      
      { if(!is.null(hlines)) geom_hline(linewidth = 0.25, color = "grey70", yintercept = hlines) } +
      
      gg_geom_ribbon +
      gg_geom_line +
      gg_geom_segment +
      
      gg_scale_color_manual +
      gg_scale_fill_manual +
      
      gg_scale_y_continuous +
      theme(legend.position = "none")
}

PLOTS[[FEATURE_TYPE$COVERAGE_DISTRIBUTION]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$COVERAGE_DISTRIBUTION)
   
   plot_labels <- labs(title = "Coverage", x = "Coverage", y = "Prop. of bases")
   
   if(is.null(plot_data)){
      return(plot_missing_data(plot_labels))
   }
   
   plot_distribution(plot_data, x = "ReadDepth", mark_sample_peak = TRUE, invert_normal = TRUE, hlines = 0) +
      plot_labels +
      render_now()
})

PLOTS[[FEATURE_TYPE$FRAG_LENGTH_DISTRIBUTION]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$FRAG_LENGTH_DISTRIBUTION)
   
   plot_labels <- labs(title = "Fragment length", x = "Fragment length", y = "Prop. of fragments")
   
   if(is.null(plot_data)){
      return(plot_missing_data(plot_labels))
   }
   
   plot_distribution(plot_data, x = "FragLength", mark_sample_peak = TRUE, invert_normal = TRUE, hlines = 0) +
      plot_labels +
      render_now()
})

PLOTS[[FEATURE_TYPE$GC_BIAS]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$GC_BIAS)
   
   plot_labels <- labs(title = "GC bias", x = "GC percentage", y = "Read depth")
   
   if(is.null(plot_data)){
      return(plot_missing_data(plot_labels))
   }
   
   plot_distribution(plot_data, x = "GCBucket", mark_sample_peak = FALSE, invert_normal = FALSE) +
      plot_labels +
      render_now()
})

## =============================
## Box plots
## =============================
PAIRWISE_PLOT_TYPE <- list(
   BOX = "box",
   POINT_RANGE = "point_range",
   BAR = "bar"
)

plot_pairwise_comparison <- function(
   plot_data, x, y = "FeatureValue", plot_type = PAIRWISE_PLOT_TYPE$BOX, hlines = NULL, vlines = NULL, box_width_scale = 1, bar_baseline = 0
){
   
   if(FALSE){
      y = "FeatureValue"
      plot_type = PAIRWISE_PLOT_TYPE$BOX
      hlines = NULL
      vlines = NULL
      box_width_scale = 1
      bar_baseline = 0
      
      x = "ReadCount"; plot_data = get_plot_data(FEATURE_TYPE$DUPLICATE_FREQ)
      x = "OriginalQualBin"; plot_data = get_plot_data(FEATURE_TYPE$BQR_BY_ORIG_QUAL); gg_facet = facet_grid("ReadType ~ StandardMutation")
   }
   
   valid_plot_types <- unlist(PAIRWISE_PLOT_TYPE, use.names = FALSE)
   if(!(plot_type %in% valid_plot_types)){
      LOGGER$error("`plot_type` must be one of: %s", paste(valid_plot_types, collapse = ", "))
   }
   
   gg_geom_point <- geom_blank()
   gg_geom_boxplot <- geom_blank()
   gg_geom_linerange <- geom_blank()
   gg_geom_bar <- geom_blank()
   gg_position_dodge <- position_dodge(width = 0.5, preserve = "single")
   
   sample_type_colors <- sapply(SAMPLE_TYPE, `[[`, "color")
   gg_scale_fill_manual <- scale_fill_manual(values = sample_type_colors, drop = FALSE)
   gg_scale_color_manual <- scale_color_manual(values = sample_type_colors, drop = FALSE)
   
   if(plot_type == PAIRWISE_PLOT_TYPE$BOX){
      gg_geom_point <- geom_point(shape = 21, position = gg_position_dodge, size = 1.8)
      
      gg_geom_boxplot <- geom_boxplot(
         aes(ymin = PctMin, lower = PctLower, middle = PctMid, upper = PctUpper, ymax = PctMax),
         position = gg_position_dodge, width = 0.4 * box_width_scale,
         stat = "identity", alpha = 0.3, size = 0.25, color = "grey70"
      )
      
      gg_scale_color_manual <- geom_blank()
   }
   
   if(plot_type == PAIRWISE_PLOT_TYPE$POINT_RANGE){ ## Used when boxplot would be too cluttered
      
      gg_geom_point <- geom_point(aes(color = SampleType), shape = 21, position = gg_position_dodge, size = 0.8)
      
      gg_geom_linerange <- geom_linerange(
         aes(color = SampleType, ymin = PctMin, ymax = PctMax),
         position = gg_position_dodge, linewidth = 0.3, linetype = "11"
      )
   }
   
   if(plot_type == PAIRWISE_PLOT_TYPE$BAR){
      
      if(!is.na(COHORT_PERCENTILES_FILE)){
         LOGGER$error("plot_type '%s' not allowed when cohort data is available", PAIRWISE_PLOT_TYPE$BAR)
      }
      
      ## Use geom_crossbar instead of geom_bar here so that bars are not inverted when values are <1 when in log scale
      gg_geom_bar <- geom_crossbar(
         aes(ymin = bar_baseline, ymax = .data[[y]]), position = gg_position_dodge,
         color = "black", middle.color = NA, linewidth = 0.3, width = 0.5*box_width_scale
      )
      
      gg_scale_color_manual <- geom_blank()
   }
   
   if(is.na(COHORT_PERCENTILES_FILE)){
      gg_geom_boxplot <- geom_blank()
      gg_geom_linerange <- geom_blank()
   }
   
   ggplot(plot_data, aes(x = .data[[x]], y = .data[[y]], fill = SampleType)) +
      
      { if(!is.null(hlines)) geom_hline(linewidth = 0.25, color = "grey70", yintercept = hlines) } +
      { if(!is.null(vlines)) geom_vline(linewidth = 0.25, color = "grey70", xintercept = vlines) } +
      
      gg_geom_boxplot + 
      gg_geom_linerange +
      gg_geom_bar +
      gg_geom_point +
      
      gg_scale_fill_manual +
      gg_scale_color_manual +
      
      # gg_facet +
      
      theme(
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         panel.spacing.x = unit(-0.5, "pt"),
         axis.text.y.right = element_blank(),
         axis.ticks.y.right = element_blank(),
         legend.position = "none"
      )
}

box_or_bar_plot <- function(){
   if(is.na(COHORT_PERCENTILES_FILE)){ 
      PAIRWISE_PLOT_TYPE$BAR
   } else { 
      PAIRWISE_PLOT_TYPE$BOX
   }
}

PLOTS[[FEATURE_TYPE$DUPLICATE_FREQ]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$DUPLICATE_FREQ)
   
   plot_labels <- labs(title = "Duplicate frequency", x = "Duplicate read count", y = "Prop. of read groups")
   
   if(is.null(plot_data)){
      return(plot_missing_data(plot_labels))
   }
   
   plot_data <- plot_data %>% dplyr::mutate(ReadCount = reverse_levels(ReadCount))
   
   plot_pairwise_comparison(plot_data, x = "ReadCount", plot_type = box_or_bar_plot()) +
      plot_labels +
      coord_flip() +
      theme(
         panel.grid.major.x = element_line(color = "grey90", linewidth = 0.25),
      ) +
      render_now()
})

PLOTS[[FEATURE_TYPE$DISCORDANT_FRAG_FREQ]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$DISCORDANT_FRAG_FREQ)
   
   plot_labels <- labs(title = "Discordant fragment frequency", x = "Discordant fragment type", y = "Prop. of reads")
   
   if(is.null(plot_data)){
      return(plot_missing_data(plot_labels))
   }
   
   plot_data <- plot_data %>% dplyr::mutate(DiscordantFragType = reverse_levels(DiscordantFragType))
   
   plot_pairwise_comparison(plot_data, x = "DisplayName", plot_type = box_or_bar_plot()) + 
      scale_y_continuous(transform = "log10", labels = function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)) +
      plot_labels +
      coord_flip() +
      theme(
         panel.grid.major.x = element_line(color = "grey90", linewidth = 0.25),
      ) +
      render_now()
})

PLOTS[[FEATURE_TYPE$MISSED_VARIANT_LIKELIHOOD]] <- local({
   
   MIN_MISSED_VARIANT_LIKELIHOOD <- 0.01
   TOP_N_GENES <- 20
   
   plot_labels <- labs(
      title = "Genes with potential missed variants",
      x = sprintf("Genes (top %s per sample type)", TOP_N_GENES), 
      y = sprintf("Missed variant likelihood (>%s)", MIN_MISSED_VARIANT_LIKELIHOOD)
   )
   
   plot_data <- get_plot_data(FEATURE_TYPE$MISSED_VARIANT_LIKELIHOOD)
   
   if(is.null(plot_data)){
      return(plot_missing_data(plot_labels))
   }
   
   plot_data <- plot_data %>% 
      dplyr::filter(FeatureValue >= MIN_MISSED_VARIANT_LIKELIHOOD) %>%
      
      dplyr::group_by(SampleType) %>% 
      dplyr::arrange(-FeatureValue) %>%
      dplyr::slice_head(n = TOP_N_GENES) %>%
      dplyr::ungroup() %>%
      
      dplyr::mutate(
         Gene = Gene %>% preordered_factor() %>% reverse_levels(),
         SampleType = SampleType %>% preordered_factor() %>% reverse_levels()
      )
   
   unique_genes <- unique(as.character(plot_data$Gene))
   gene_count <- length(unique_genes)

   ## Split plot into 2 panels if there are too many genes
   N_GENES_TO_SPLIT_PANEL <- 15
   gg_facet_wrap <- geom_blank()
   if(gene_count > N_GENES_TO_SPLIT_PANEL){
      midpoint <- round(gene_count / 2)
      gene_group <- as.integer(1:gene_count > midpoint); names(gene_group) <- unique_genes
      
      plot_data$GeneGroup <- gene_group[as.character(plot_data$Gene)]
      gg_facet_wrap <- facet_wrap(". ~ GeneGroup", scales = "free_y")
   }
   
   plot_pairwise_comparison(plot_data, x = "Gene", plot_type = box_or_bar_plot()) + 
      plot_labels +
      gg_facet_wrap +
      scale_y_continuous(labels = scales::label_number(drop0trailing=TRUE)) +
      coord_flip(ylim = c(0, NA)) +
      theme(
         panel.grid.major.x = element_line(color = "grey90", linewidth = 0.25),
         plot.title = element_text(hjust = 0.5),
         axis.title.y = element_text(size = 10),
         axis.text.y = element_text(size = 8),
         panel.spacing.x = unit(3, "pt"),
         strip.text.x = element_blank(),
         strip.background = element_blank(),
      ) +
      render_now()
})

PLOTS[[FEATURE_TYPE$BQR_BY_ORIG_QUAL]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$BQR_BY_ORIG_QUAL)
   
   plot_labels <- labs(title = "BQR by original base quality", x = "Original base quality", y = "Phred score adjustment")
   
   if(is.null(plot_data)){
      return(plot_missing_data(plot_labels))
   }
   
   plot_pairwise_comparison(plot_data, x = "OriginalQualBin", hlines = 0) +
      facet_grid("ReadType ~ StandardMutation") +
      scale_y_continuous(
         labels = function(x){ ifelse(x > 0, paste0("+",x), x) },
         sec.axis = dup_axis(name = "Consensus type")
      ) +
      plot_labels + 
      theme(
         axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 7),
      ) +
      render_now()
})

PLOTS[[FEATURE_TYPE$BQR_BY_SNV96_CONTEXT]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$BQR_BY_SNV96_CONTEXT)
   
   plot_labels <- labs(title = "BQR by SNV96 context", x = "Mutation context", y = "Phred score adjustment")
   
   if(is.null(plot_data)){
      plot_labels$title <- paste0(plot_labels$title, ", base quality: ", plot_data$OriginalQualBin[1])
      return(plot_missing_data(plot_labels))
   }
   
   plot_pairwise_comparison(
      plot_data, x = "StandardTrinucContext", plot_type = PAIRWISE_PLOT_TYPE$POINT_RANGE, 
      hlines = 0, vlines = c(4, 8, 12) + 0.5
   ) +
      facet_grid("ReadType ~ StandardMutation", scales = "free_x") +
      scale_y_continuous(
         labels = function(x){ ifelse(x > 0, paste0("+",x), x) },
         sec.axis = dup_axis(name = "Consensus type")
      ) +
      plot_labels +
      theme(
         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      ) +
      render_now()
})

PLOTS[[FEATURE_TYPE$MS_INDEL_ERROR_RATES]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$MS_INDEL_ERROR_RATES)
   
   plot_labels <- labs(title = "Microsatellite indel error rates", x = "Repeat units", y = "Phred score") 
   
   if(is.null(plot_data)){
      return(plot_missing_data(plot_labels))
   }
   
   plot_data$RefNumUnits <- as.numeric(plot_data$RefNumUnits)

   plot_pairwise_comparison(plot_data, x = "RefNumUnits", plot_type = PAIRWISE_PLOT_TYPE$POINT_RANGE) +
      facet_grid("ConsensusType ~ RepeatUnitType") +
      scale_x_continuous(
         breaks = scales::breaks_width(3), 
         limits = range(c(3, 12, plot_data$RefNumUnits), na.rm = TRUE)
      ) +
      scale_y_continuous(limits = c(0, NA), sec.axis = dup_axis(name = "Consensus type")) +
      plot_labels +
      theme(
         panel.grid.major = element_line(color = "grey90", linewidth = 0.25)
      ) +
      render_now()
})

PLOTS[[FEATURE_TYPE$MS_INDEL_ERROR_BIAS]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$MS_INDEL_ERROR_BIAS)
   
   plot_labels <- labs(
      title = "Microsatellite indel error bias",
      x = "Repeat units",
      y = expression(atop("Phred score diff.",atop("more del. errors <-> more ins. errors")))
   )
   
   if(is.null(plot_data)){
      return(plot_missing_data(plot_labels))
   }
   
   plot_data$RefNumUnits <- as.numeric(plot_data$RefNumUnits)
   
   plot_pairwise_comparison(plot_data, x = "RefNumUnits", plot_type = PAIRWISE_PLOT_TYPE$POINT_RANGE) +
      facet_grid("ConsensusType ~ RepeatUnitType") +
      scale_x_continuous(
         breaks = scales::breaks_width(3), 
         limits = range(c(3, 12, plot_data$RefNumUnits), na.rm = TRUE)
      ) +
      scale_y_continuous(
         labels = function(x){ ifelse(x > 0, paste0("+",x), x) },
         sec.axis = dup_axis(name = "Consensus type")
      ) +
      plot_labels +
      theme(
         panel.grid.major = element_line(color = "grey90", linewidth = 0.25)
      ) +
      render_now()
})

## =============================
## Summary table
## =============================

FEATURE_GROUP <- list(
   MAPPING = "Mapping",
   COPY_NUMBER = "Copy number",
   CONTAMINATION = "Contamination",
   MUTATIONAL_BURDEN = "Mutational burden"
)

NUMBER_FORMAT <- list(
   NUMBER = "NUMBER",
   PERCENT = "PERCENT",
   LOG10 = "LOG10"
)

SUMMARY_TABLE_DATA <- get_plot_data(FEATURE_TYPE$SUMMARY_TABLE)

get_sub_table_data <- function(feature_group, number_format){

   plot_data <- SUMMARY_TABLE_DATA %>%
      dplyr::filter(FeatureGroup == feature_group & NumberFormat == number_format) %>%
      dplyr::mutate(
         DisplayName = reverse_levels(DisplayName),
         SampleType = reverse_levels(SampleType)
      )

   attr(plot_data, "number_format") <- number_format
   attr(plot_data, "feature_group") <- feature_group

   return(plot_data)
}

get_limits <- function(plot_data, minimal_limits = c(NA, NA)){
   
   #' Ensure y axis limits are wide enough to show pretty breaks.
   #' 
   #' Given:
   #'   minimal limits : A=====B
   #'   range in cohort:    C====D
   #'   value in sample:          E
   #'
   #' Chosen limits are: A and E
   #'
   #' For example, the min/max TMB (based on  sample and/or cohort data) could be 0.5 to 1, but we
   #' know TMB can go form 0 to 100s/1000s. We would therefore so the minimal limits to
   #' e.g. c(0, 1000).
   #'
   #' Providing NA minimal limit values lets ggplot decide the limits based on the min/max of the data
   #' 
   
   values <- c(minimal_limits, plot_data$FeatureValue, plot_data$PctMin, plot_data$PctMax)
   limits <- range(values, na.rm = TRUE)
   
   ## Let ggplot handle limits with infinite values
   limits[!is.finite(limits)] <- NA 
   
   ## User intends to let ggplot handle limits
   if(is.na(minimal_limits[1]))
      limits[1] <- NA
   
   if(is.na(minimal_limits[2]))
      limits[2] <- NA
   
   return(limits)
}

get_div_lines <- function(n_table_rows, direction = "horizontal"){
   
   positions <- seq(1.5, n_table_rows)
   color <- if(n_table_rows > 1) "grey70" else "white"
   linetype <- "dotted"
   
   if(direction == "horizontal"){
      geom_hline(yintercept = positions, color = color, linetype = linetype)
   } else if(direction == "vertical"){
      geom_vline(xintercept = positions, color = color, linetype = linetype)
   } else {
      LOGGER$error("`direction` must be 'horizontal' or 'vertical'")
   }
}

get_qc_status_enums <- function(qc_status_strings){
   #qc_status_strings = SUMMARY_TABLE_DATA$QcStatus
   qc_status_strings <- as.character(qc_status_strings)
   
   qc_status_enums <- dplyr::case_when(
      startsWith(qc_status_strings, QC_STATUS$FAIL$name) ~ QC_STATUS$FAIL$name,
      startsWith(qc_status_strings, QC_STATUS$WARN$name) ~ QC_STATUS$WARN$name,
      .default = QC_STATUS$PASS$name
   )
   
   factor(qc_status_enums, sapply(QC_STATUS, `[[`, "name"))
}

plot_sub_table <- function(plot_data, show_title = FALSE, show_sample_type_label = FALSE, axis_limits = waiver()){

   if(FALSE){
      show_title = TRUE
      show_sample_type_label = TRUE
      axis_limits = c(NA, NA)
      
      plot_data = get_sub_table_data(FEATURE_GROUP$MAPPING, NUMBER_FORMAT$PERCENT)
      axis_limits = c(0,1)
      
      plot_data = get_sub_table_data(FEATURE_GROUP$MAPPING, NUMBER_FORMAT$NUMBER)
      axis_limits = get_limits(plot_data, c(0, 100))
      
      plot_data = get_sub_table_data(FEATURE_GROUP$COPY_NUMBER, NUMBER_FORMAT$LOG10)
      axis_limits = get_limits(plot_data, c(NA, 1000))
      
      plot_data = get_sub_table_data(FEATURE_GROUP$MUTATIONAL_BURDEN, NUMBER_FORMAT$LOG10)
      axis_limits = get_limits(plot_data, c(NA, 1000))
   }

   feature_group <- attr(plot_data, "feature_group")
   number_format <- attr(plot_data, "number_format")
   n_rows <- plot_data$DisplayName %>% unique() %>% length()

   LOGGER$debug("  featureGroup(%s) numberFormat(%s)", feature_group, number_format)
      
   ## Feature values =============================
   if(number_format == NUMBER_FORMAT$PERCENT){
      value_fmt_func <- scales::label_percent(accuracy = 0.1)
   } else {
      value_fmt_func <- function(x) ifelse(x < 1, signif(x, 2), round(x, 1))
   }
   
   ## Prep data
   plot_data <- plot_data %>% dplyr::mutate(
      QcStatusEnum = get_qc_status_enums(QcStatus),
      QcStatusEnum = factor(QcStatusEnum, sapply(QC_STATUS, `[[`, "name"))
   )
   
   ## Aesthetics
   qc_status_colors <- sapply(QC_STATUS, `[[`, "color")
   sample_type_colors <- sapply(SAMPLE_TYPE, `[[`, "color")
   sample_type_names <- sapply(SAMPLE_TYPE, `[[`, "display_name")
   
   ## Plot
   subplot_values <- ggplot(plot_data, aes(y = DisplayName, x = SampleType, group = SampleType, color = SampleType)) +
      
      geom_label(
         aes(label = value_fmt_func(FeatureValue), fill = QcStatusEnum), 
         size = 3, label.padding = unit(4, "pt"),
         border.colour = ifelse(plot_data$QcStatusEnum == QC_STATUS$PASS$name, "#FFFFFF00", "black")
      ) +
      scale_color_manual(values = sample_type_colors) +
      scale_fill_manual(values = qc_status_colors) +
      
      get_div_lines(n_rows, "horizontal") +
      
      scale_x_discrete(position = "bottom", limits = rev, drop = FALSE, labels = sample_type_names) +
      guides(color = "none") +
      labs(title = feature_group) +
      
      theme(
         plot.title = if(show_title) element_text() else element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.ticks.x = if(show_sample_type_label) element_line() else element_blank(),
         axis.text.x = if(show_sample_type_label) element_text() else element_blank(),
         plot.margin = margin(b = 10),
         legend.position = "none"
      )

   ## Sample vs cohort =============================
   box_width_scale <- length(unique(plot_data$SampleType)) / length(levels(plot_data$SampleType))
   plot_type <- if(is.na(COHORT_PERCENTILES_FILE)) PAIRWISE_PLOT_TYPE$BAR else PAIRWISE_PLOT_TYPE$BOX
   
   subplot_pairwise_comparison <- 
      plot_pairwise_comparison(plot_data, x = "DisplayName", y = "FeatureValue", box_width_scale = box_width_scale, plot_type = plot_type) +
      get_div_lines(n_rows, "vertical") + 
      scale_y_continuous(
         
         limits = axis_limits,
         transform = if(number_format == NUMBER_FORMAT$LOG10) "log10" else "identity",
         
         label = if(number_format == NUMBER_FORMAT$PERCENT){ 
            scales::label_percent()
         } else if(number_format == NUMBER_FORMAT$LOG10) {
            function(x) format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE)
         } else {
            waiver()
         }
      ) +
      coord_flip() +
      theme(
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.line.y = element_blank(),
         plot.margin = margin(b = 10, r = 15, l = 2),
         legend.position = "none"
      ) +
      render_now()

   ## Combine plots =============================
   subplots_combined <- patchwork::wrap_plots(subplot_values, subplot_pairwise_comparison, nrow = 1)
   
   subplots_combined$height <- n_rows
   
   subplots_combined
}

PLOTS[[FEATURE_TYPE$SUMMARY_TABLE]] <- local({
   
   plots <- list()
   
   ## Mapping ================================
   plots[[1]] <- 
      get_sub_table_data(FEATURE_GROUP$MAPPING, NUMBER_FORMAT$NUMBER) %>% 
      plot_sub_table(show_title = TRUE, axis_limits = get_limits(., minimal_limits = c(0, 100)))
   
   plots[[2]] <- 
      get_sub_table_data(FEATURE_GROUP$MAPPING, NUMBER_FORMAT$PERCENT) %>% 
      plot_sub_table(axis_limits = c(0, 1))
   
   ## Copy number ================================
   plots[[3]] <- 
      get_sub_table_data(FEATURE_GROUP$COPY_NUMBER, NUMBER_FORMAT$NUMBER) %>% 
      plot_sub_table(show_title = TRUE, axis_limits = c(0, NA))
   
   plots[[4]] <- 
      get_sub_table_data(FEATURE_GROUP$COPY_NUMBER, NUMBER_FORMAT$PERCENT) %>% 
      plot_sub_table(axis_limits = c(0, 1))
   
   plots[[5]] <- 
      get_sub_table_data(FEATURE_GROUP$COPY_NUMBER, NUMBER_FORMAT$LOG10) %>% 
      plot_sub_table(axis_limits = get_limits(., minimal_limits =c(NA, 1000)))
   
   ## Contamination ================================
   plots[[6]] <- get_sub_table_data(FEATURE_GROUP$CONTAMINATION, NUMBER_FORMAT$PERCENT) %>% 
      plot_sub_table(show_title = TRUE, axis_limits = c(0, 1))
   
   ## Mutational burden ================================
   plots[[7]] <- 
      get_sub_table_data(FEATURE_GROUP$MUTATIONAL_BURDEN, NUMBER_FORMAT$LOG10) %>% 
      plot_sub_table(show_title = TRUE, show_sample_type_label = TRUE, axis_limits = get_limits(., minimal_limits = c(NA, 1000)))
   
   heights <- sapply(plots, function(p){ p$height })
   
   plots_combined <- patchwork::wrap_plots(plots, ncol = 1, heights = heights)
   plots_combined
})

## =============================
## Other plot components
## =============================

PLOT_NAME_LEGEND <- "CUSTOM_LEGEND"

PLOTS[[PLOT_NAME_LEGEND]] <- local({
   
   #' A custom legend is created for a few reasons:
   #' 
   #' 1) `patchwork::wrap_plots(..., guides = "collect")` is unreliable, especially if some plots 
   #'    for example have 'TUMOR' and 'NORMAL' sample types but other plots only have 'TUMOR' 
   #'    (but not 'NORMAL'). 
   #' 
   #' 2) We don't need to for example have the boxplot legend as well as the line plot legend both 
   #'    showing the same colors. Therefore, the custom legend only shows colors.
   #'    
   #' 3) The custom legend is also itself a plot, which allows it to be positioned in an empty spot 
   #'    with patchwork.
   
   plot_data_sample_type <- data.frame(
      SampleType = c(
         paste(SAMPLE_TYPE$TUMOR$display_name, "sample"),
         paste(SAMPLE_TYPE$TUMOR$display_name, "cohort"),
         paste(SAMPLE_TYPE$NORMAL$display_name, "sample"),
         paste(SAMPLE_TYPE$NORMAL$display_name, "cohort")
      ),
      
      Color = c(
         SAMPLE_TYPE$TUMOR$color,
         SAMPLE_TYPE$TUMOR$color %>% adjustcolor(alpha.f = 0.2),
         SAMPLE_TYPE$NORMAL$color,
         SAMPLE_TYPE$NORMAL$color %>% adjustcolor(alpha.f = 0.2)
      )
   ) %>% 
      dplyr::mutate(Index = 1:dplyr::n())
   
   plot_data_qc_status <- data.frame(
      Index = 1:length(QC_STATUS),
      QcStatus = lapply(QC_STATUS, `[[`, "name") %>% unlist(use.names = FALSE),
      Color = lapply(QC_STATUS, `[[`, "color") %>% unlist(use.names = FALSE)
   )
   
   dummy_plot <- ggplot(data.frame(), aes(SampleType, QcStatus)) +
      
      ## Sample type
      geom_point(
         data = plot_data_sample_type, 
         mapping = aes(x = Index, y = 1, color = preordered_factor(SampleType)),
         shape = 15, size = 6
      ) +
      scale_color_manual(
         name = "Sample type",
         values = plot_data_sample_type %>% dplyr::pull(Color, name = SampleType),
         labels = plot_data_sample_type %>% dplyr::pull(SampleType, name = SampleType),
      ) +
      
      ## QC status
      geom_label(
         data = plot_data_qc_status, 
         mapping = aes(x = Index, y = 2, fill = preordered_factor(QcStatus), label = QcStatus),
         color = "#FFFFFF00", border.color = "black"
      ) +
      scale_fill_manual(
         name = "QC status",
         values = plot_data_qc_status %>% dplyr::pull(Color, name = QcStatus)
      ) +
      
      ## Set legend order
      guides(
         color = guide_legend(order = 1),
         fill = guide_legend(order = 2)
      ) +
      
      theme(
         legend.position="bottom", 
         legend.direction="vertical",
      )
   
   legend <- local({ 
      plot_as_gtable <- dummy_plot %>% ggplot_build() %>% ggplot_gtable() 
      
      grobs <- plot_as_gtable$grobs
      grob_names <- sapply(grobs, function(x) x$name)
      legend_index <- which(grob_names == "guide-box")
      
      grobs[[legend_index]] %>% patchwork::wrap_elements()
   })
   
   return(legend)
})

REPORT_TITLE <- local({
   
   qc_status_enums <- get_qc_status_enums(SUMMARY_TABLE_DATA$QcStatus)
   
   form_qc_string <- function(qc_status_level, sample_type){
      df <- SUMMARY_TABLE_DATA %>% filter(
         SampleType == sample_type$name & 
         qc_status_enums == qc_status_level$name
      )
      
      if(nrow(df) == 0)
         return(character())
      
      paste0(
         qc_status_level$name, ": ",
         sprintf("%s (%s)", df$DisplayName, df$QcThreshold) %>% paste(collapse = ", ")
      )
   }
   
   form_sample_type_qc_string <- function(sample_type){
      qc_strings <- c(
         form_qc_string(QC_STATUS$FAIL, sample_type),
         form_qc_string(QC_STATUS$WARN, sample_type)
      )
      
      if(length(qc_strings) == 0){
         return(character())
      }
      
      paste(c(sample_type$display_name, qc_strings), collapse = " - ")
   }
   
   qc_strings <- c(
      form_sample_type_qc_string(SAMPLE_TYPE$TUMOR),
      form_sample_type_qc_string(SAMPLE_TYPE$NORMAL)
   )
   subtitle <- if(length(qc_strings) > 0) paste(qc_strings, collapse = "\n") else waiver()
   
   most_severe_qc_status <- qc_status_enums %>% sort() %>% head(1)
   title <- paste0(TUMOR_ID, ": ", most_severe_qc_status)
   
   patchwork::plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
   )
})

## =============================
## Combine plots
## =============================

create_report <- local({
   
   LOGGER$info("Combining plots")
   
   plot_order <- c(
      "A" = FEATURE_TYPE$SUMMARY_TABLE,
      
      "B" = FEATURE_TYPE$COVERAGE_DISTRIBUTION,
      "C" = FEATURE_TYPE$FRAG_LENGTH_DISTRIBUTION,
      "D" = FEATURE_TYPE$GC_BIAS,
      "E" = FEATURE_TYPE$DUPLICATE_FREQ,
      "F" = FEATURE_TYPE$DISCORDANT_FRAG_FREQ,
      "G" = FEATURE_TYPE$MISSED_VARIANT_LIKELIHOOD,
      
      "H" = FEATURE_TYPE$BQR_BY_ORIG_QUAL,
      "I" = FEATURE_TYPE$BQR_BY_SNV96_CONTEXT,
      "J" = FEATURE_TYPE$MS_INDEL_ERROR_RATES,
      "K" = FEATURE_TYPE$MS_INDEL_ERROR_BIAS,
      
      "Z" = PLOT_NAME_LEGEND
   )

   plots <- lapply(plot_order, function(feature_type){

      p <- PLOTS[[feature_type]]

      has_wide_y_axis <- feature_type %in% c(FEATURE_TYPE$DISCORDANT_FRAG_FREQ, FEATURE_TYPE$MISSED_VARIANT_LIKELIHOOD)
      if(has_wide_y_axis){
         p <- patchwork::free(p, "panel", side = "l")
      } else {
         p <- patchwork::free(p, "label")
      }

      return(p)
   })

   design <- "
      AABBCCDD
      AAEEFFGG
      AAHHIIII
      AAJJKKZZ
   "
   
   plots_combined <- patchwork::wrap_plots(plots, design = design) + REPORT_TITLE

   plot_width <- 20
   plot_height <- 12

   ## Suppress warnings to prevent ggplot build warnings showing again - these are already shown when calling render_now()

   LOGGER$info("Writing PDF: %s", PLOT_PATH)
   suppressWarnings({
      ggsave(filename = PLOT_PATH, plot = plots_combined, width = plot_width, height = plot_height, units = "in")
   })

   if(SINGLE_PATIENT_MODE)
   {
      ## PNGs are easier to embed in downstream reports (e.g. on a website)
      png_path <- sub("pdf$", "png", PLOT_PATH)
      LOGGER$info("Writing PNG: %s", png_path)
      suppressWarnings({
         ggsave(filename = png_path, plot = plots_combined, width = plot_width, height = plot_height, units = "in")
      })
   }
})



