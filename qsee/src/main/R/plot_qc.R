suppressPackageStartupMessages(library(dplyr))

library(ggplot2)
theme_set(
   theme_bw() +
   theme(
      panel.grid = element_blank()
   )
)

library(patchwork)
library(scales)

## =============================
## Config
## =============================

args <- commandArgs(trailingOnly = TRUE)

TUMOR_ID <- args[1]
NORMAL_ID <- ifelse(args[2] == "NA", NA, args[2])
SAMPLE_FEATURES_FILE <- args[3]
COHORT_PERCENTILES_FILE <- args[4]
OUTPUT_PATH <- args[5]
GLOBAL_LOG_LEVEL <- args[6]

if(FALSE){
    TUMOR_ID <- "TUMOR"
    NORMAL_ID <- "TUMOR-ref"

    COHORT_PERCENTILES_FILE <- "COHORT.qsee.percentiles.tsv.gz"

    output_dir <- ""
    SAMPLE_FEATURES_FILE <- sprintf("%s/%s.qsee.vis.features.tsv.gz", output_dir, TUMOR_ID)
    OUTPUT_PATH <- sprintf("%s/%s.qsee.vis.report.pdf", output_dir, TUMOR_ID)

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

   log_message <- sprintf("%s [R] [%-5s] %s", current_time, log_level$name, sprintf(fmt, ...))

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

LOGGER$debug("Running script with args:")
LOGGER$debug(" tumor_id: %s", TUMOR_ID)
LOGGER$debug(" normal_id: %s", NORMAL_ID)
LOGGER$debug(" sample_features_file: %s", SAMPLE_FEATURES_FILE)
LOGGER$debug(" cohort_percentiles_file: %s", COHORT_PERCENTILES_FILE)
LOGGER$debug(" output_path: %s", OUTPUT_PATH)
LOGGER$debug(" log_level: %s", GLOBAL_LOG_LEVEL)

## =============================
## String functions
## =============================

preordered_factor <- function(x){ factor(x, unique(x)) }

preordered_factors <- function(df){ 
   as.data.frame(lapply(df, function(x){
      if(is.character(x)){
         x <- factor(x, unique(x))
      }
      return(x)
   }))
}

KEY_VALUE_SEP <- "="
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

## =============================
## Constants
## =============================

SAMPLE_TYPE <- list(
   TUMOR = list(name = "TUMOR", human_readable_name = "Tumor", color = "#D1392C"),
   NORMAL = list(name = "NORMAL", human_readable_name = "Normal", color = "#4A7DB4")
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
   LIM_LOWER = list(name = "PctLimLower", colname = paste0(PERCENTILE_PREFIX, 0)),
   
   MIN   = list(name = "PctMin"  , colname = paste0(PERCENTILE_PREFIX,  5)),
   LOWER = list(name = "PctLower", colname = paste0(PERCENTILE_PREFIX, 25)),
   MID   = list(name = "PctMid"  , colname = paste0(PERCENTILE_PREFIX, 50)),
   UPPER = list(name = "PctUpper", colname = paste0(PERCENTILE_PREFIX, 75)),
   MAX   = list(name = "PctMax"  , colname = paste0(PERCENTILE_PREFIX, 95)),
   
   LIM_UPPER = list(name = "PctLimUpper", colname = paste0(PERCENTILE_PREFIX, 100))
)

NUMBER_FORMATS <- list(
   NUMBER = "NUMBER",
   PERCENT = "PERCENT",
   LOG = "LOG"
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
   LOGGER$info("Loading sample features from: %s", SAMPLE_FEATURES_FILE)
   
   sample_features <- read.delim(SAMPLE_FEATURES_FILE, na.strings = c("NA", "null"))
   sample_features <- sample_features %>% dplyr::filter(SampleId %in% c(TUMOR_ID, NORMAL_ID))
   
   return(sample_features)
}

COHORT_DATA <- load_cohort_percentiles()
SAMPLE_DATA <- load_sample_features()

PLOTS <- list()

## =============================
## Plot helper functions
## =============================

plot_missing_data <- function(plot_labels = labs()){

   ggplot() +
      annotate("text", x = 0, y = 0, label = "Missing sample data") +
      plot_labels +
      theme(
         axis.text = element_blank(),
         axis.ticks = element_blank()
      )
}

get_plot_data <- function(feature_type){
   
   LOGGER$info("Plotting featureType(%s)", feature_type)
   
   ## Select rows
   cohort_data <- COHORT_DATA %>% dplyr::filter(FeatureType == feature_type)
   sample_data <- SAMPLE_DATA %>% dplyr::filter(FeatureType == feature_type)
   
   if(nrow(sample_data) == 0){
      return(data.frame())
   }
   
   merged_data <- merge(
      sample_data, cohort_data, 
      by=c("SampleType", "SourceTool", "FeatureType", "FeatureName"), 
      all.x = TRUE, ## Only select the features in the cohort data that are present in the sample
      sort = FALSE
   )
   
   ## Split multi-field strings
   for(colname in colnames(merged_data)){
      column <- merged_data[[colname]]
      
      if(!has_multiple_fields(column))
         next
      
      column_as_df <- strings_to_df(column)
      
      merged_data[[colname]] <- NULL
      merged_data <- cbind(merged_data, column_as_df)
   }

   merged_data <- preordered_factors(merged_data)
   
   return(merged_data)
}

## =============================
## Box plots
## =============================

plot_boxplot <- function(
   plot_data, x, y = "FeatureValue", plot_labels = geom_blank(), plot_type = "boxplot",
   hlines = NULL, vlines = NULL
){
   
   if(FALSE){
      plot_type = "boxplot"
      hlines = NULL
      vlines = NULL
      
      x = "OriginalQualBin"
      y = "FeatureValue"

      plot_data <- get_plot_data(FEATURE_TYPE$BQR_BY_ORIG_QUAL$name)
      plot_labels <- labs(title = "BQR by original base quality", x = "Original base quality", y = "Phred score adjustment")
      gg_facet <- facet_grid("ReadType ~ StandardMutation")
   }
   
   if(nrow(plot_data) == 0){
      return(plot_missing_data(plot_labels))
   }
   
   sample_type_colors <- sapply(SAMPLE_TYPE, `[[`, "color")
   gg_position_dodge <- position_dodge(width = 0.5)
   boxplot_width <- 0.2 * (plot_data$SampleType %>% unique() %>% length())
   
   geom_sample <- geom_blank()
   geom_cohort <- geom_blank()
   
   if(plot_type == "boxplot"){
      geom_sample <- geom_point(
         shape = 21, position = gg_position_dodge, size = 1.5
      )
      
      geom_cohort <- geom_boxplot(
         aes(ymin = PctMin, lower = PctLower, middle = PctMid, upper = PctUpper, ymax = PctMax),
         position = gg_position_dodge, width = boxplot_width,
         stat = "identity", alpha = 0.3, size = 0.25, color = "grey70"
      )
   } else if(plot_type == "pointrange") {
      geom_sample <- geom_point(
         aes(color = SampleType), 
         shape = 21, position = gg_position_dodge, size = 1
      )
      
      geom_cohort <- geom_linerange(
         aes(color = SampleType, ymin = PctMin, ymax = PctMax),
         position = gg_position_dodge, linewidth = 0.3, linetype = "11"
      )
   } else {
      LOGGER$error("`plot_type` must be 'boxplot' or 'pointrange'")
   }
   
   ggplot(plot_data, aes(x = .data[[x]], y = .data[[y]], fill = SampleType)) + 
      
      { if(!is.null(hlines)) geom_hline(linewidth = 0.25, color = "grey70", yintercept = hlines) } +
      { if(!is.null(vlines)) geom_vline(linewidth = 0.25, color = "grey70", xintercept = vlines) } +
      
      geom_cohort +
      geom_sample +
      
      scale_fill_manual(values = sample_type_colors) +
      scale_color_manual(values = sample_type_colors) +
      
      #gg_facet +
      
      plot_labels +
      
      theme(
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         panel.spacing.x = unit(-0.5, "pt"),
         axis.text.y.right = element_blank(),
         axis.ticks.y.right = element_blank(),
         legend.position = "none"
      )
}

PLOTS[[FEATURE_TYPE$DUPLICATE_FREQ]] <- local({
   plot_data <- get_plot_data(FEATURE_TYPE$DUPLICATE_FREQ)
   
   plot_labels <- labs(title = "Duplicate frequency", x = "Duplicate read count", y = "Prop. of read groups")
   
   plot_boxplot(plot_data, x = "ReadCount", plot_labels = plot_labels) +
      theme(
         panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
         axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
      )
})

PLOTS[[FEATURE_TYPE$DISCORDANT_FRAG_FREQ]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$DISCORDANT_FRAG_FREQ)
   
   plot_labels <- labs(
      title = "Discordant fragment frequency", 
      x = "Discordant fragment type", 
      y = "Prop. of reads"
   )
   
   plot_boxplot(plot_data, x = "DiscordantFragType", plot_labels = plot_labels) + 
      scale_y_continuous(
         transform = "log10", 
         labels = scales::trans_format("log10", scales::math_format(10^.x))
      ) +
      theme(
         axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      )
})

PLOTS[[FEATURE_TYPE$MISSED_VARIANT_LIKELIHOOD]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$MISSED_VARIANT_LIKELIHOOD)
   
   MIN_MISSED_VARIANT_LIKELIHOOD <- 0.01
   TOP_N_GENES <- 20
   
   plot_labels <- labs(
      title = sprintf("Top %s genes with potential missed variants", TOP_N_GENES),
      x = "Gene", 
      y = "Missed variant likelihood"
   )
   
   if(nrow(plot_data) == 0){
      return(plot_missing_data(plot_labels))
   }
   
   sample_genes_of_interest <- plot_data %>% 
      dplyr::filter(
         FeatureValue >= MIN_MISSED_VARIANT_LIKELIHOOD | 
         PctMid >= MIN_MISSED_VARIANT_LIKELIHOOD
      ) %>%
      dplyr::arrange(-FeatureValue) %>%
      dplyr::pull(Gene) %>% 
      unique() %>% 
      head(TOP_N_GENES)
   
   plot_data <- plot_data %>% 
      dplyr::filter(Gene %in% sample_genes_of_interest) %>% 
      dplyr::mutate(Gene = factor(Gene, sample_genes_of_interest))
   
   plot_boxplot(plot_data, x = "Gene", plot_labels = plot_labels) + 
      theme(
         axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      )
})

PLOTS[[FEATURE_TYPE$BQR_BY_ORIG_QUAL]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$BQR_BY_ORIG_QUAL)
   
   plot_labels <- labs(
      title = "BQR by original base quality", 
      x = "Original base quality", 
      y = "Phred score adjustment"
   )
   
   plot_boxplot(plot_data, x = "OriginalQualBin", plot_labels = plot_labels, hlines = 0) +
      facet_grid("ReadType ~ StandardMutation") +
      scale_y_continuous(
         labels = function(x){ ifelse(x > 0, paste0("+",x), x) },
         sec.axis = dup_axis(name = "Consensus type")
      ) +
      theme(
         axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      )
})

PLOTS[[FEATURE_TYPE$BQR_BY_SNV96_CONTEXT]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$BQR_BY_SNV96_CONTEXT)
   
   plot_labels <- labs(title = "BQR by SNV96 context", x = "Mutation context", y = "Phred score adjustment")
   
   if(nrow(plot_data) > 0){
      plot_labels$title <- paste0(plot_labels$title, ", base quality: ", plot_data$OriginalQualBin[1])
   }
   
   plot_boxplot(
      plot_data, x = "StandardTrinucContext", plot_labels = plot_labels, plot_type = "pointrange", 
      hlines = 0, vlines = c(4, 8, 12) + 0.5
   ) +
      facet_grid("ReadType ~ StandardMutation", scales = "free_x") +
      scale_y_continuous(
         labels = function(x){ ifelse(x > 0, paste0("+",x), x) },
         sec.axis = dup_axis(name = "Consensus type")
      ) +
      theme(
         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
      )
})

PLOTS[[FEATURE_TYPE$MS_INDEL_ERROR_RATES]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$MS_INDEL_ERROR_RATES)
   
   plot_labels <- labs(title = "Microsatellite indel error rates", x = "Repeat units", y = "Phred score") 
   
   plot_boxplot(plot_data, plot_labels = plot_labels, x = "RefNumUnits", plot_type = "pointrange") +
      facet_grid("ConsensusType ~ RepeatUnitType") +
      scale_x_discrete(breaks = function(x) ifelse(as.numeric(x) %% 3 == 0, x, "") ) +
      scale_y_continuous(limits = c(0, NA), sec.axis = dup_axis(name = "Consensus type")) +
      theme(
         panel.grid.major = element_line(color = "grey90", linewidth = 0.25)
      )
})

PLOTS[[FEATURE_TYPE$MS_INDEL_ERROR_BIAS]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$MS_INDEL_ERROR_BIAS)
   
   plot_labels <- labs(title = "Microsatellite indel error bias", x = "Repeat units", y = "Phred score diff.")
   
   if(nrow(plot_data) == 0){
      return(plot_missing_data(plot_labels))
   }
   
   text_data <- data.frame(
      SampleType     = SAMPLE_TYPE$TUMOR$name,
      RepeatUnitType = plot_data$RepeatUnitType[1],
      RefNumUnits    = c(-Inf, -Inf),
      FeatureValue   = c( Inf, -Inf),
      Label          = c("More ins. errors", "More del. errors"),
      Vjust          = c(1.5, -0.5)
   )
   
   plot_boxplot(plot_data, x = "RefNumUnits", plot_labels = plot_labels, hlines = 0, plot_type = "pointrange") +
      facet_grid("ConsensusType ~ RepeatUnitType") +
      geom_text(
         data = text_data, mapping = aes(label = Label, vjust = Vjust), 
         hjust = -0.05, size = 2.5
      ) +
      scale_x_discrete(breaks = function(x) ifelse(as.numeric(x) %% 3 == 0, x, "") ) +
      scale_y_continuous(
         labels = function(x){ ifelse(x > 0, paste0("+",x), x) },
         sec.axis = dup_axis(name = "Consensus type")
      )
})

## =============================
## Line / PDF
## =============================

plot_distribution <- function(plot_data, x, plot_labels = geom_blank(), invert_normal = FALSE, mark_sample_peak = FALSE){
   
   if(FALSE){
      plot_data = get_plot_data(FEATURE_TYPE$COVERAGE_DISTRIBUTION)
      plot_labels = labs(title = "Coverage", x = "Coverage", y = "Prop. of bases")
      x = "ReadDepth"
   }
   
   if(nrow(plot_data) == 0){
      return(plot_missing_data(plot_labels))
   }
   
   plot_data[[x]] <- as.numeric(plot_data[[x]])
   
   gg_hline <- geom_blank()
   gg_scale_y_continuous <- geom_blank()
   if(invert_normal){
      
      signs <- ifelse(plot_data$SampleType == SAMPLE_TYPE$NORMAL$name, -1, 1)
      
      plot_data <- plot_data %>% mutate(
         PctMin   = signs * PctMin,
         PctLower = signs * PctLower,
         PctMid   = signs * PctMid,
         PctUpper = signs * PctUpper,
         PctMax   = signs * PctMax,
         FeatureValue = signs * FeatureValue
      )
      
      gg_hline <- geom_hline(yintercept = 0, linewidth = 0.25)
      gg_scale_y_continuous <- scale_y_continuous(labels = function(x) abs(x))
   }
   
   geom_sample_peak <- geom_blank()
   if(mark_sample_peak){
      
      peak_data <- plot_data %>% 
         dplyr::group_by(SampleType) %>%
         
         ## Remove the bottom 5% since the coverage plot has a spike there
         dplyr::slice(floor(dplyr::n() * 0.05):ceiling(dplyr::n() * 0.95)) %>% 
         
         dplyr::reframe(
            XPos = .data[[x]][which.max(abs(FeatureValue))],
            Height = FeatureValue[which.max(abs(FeatureValue))]
         )
      
      geom_sample_peak <- geom_segment(
         data = peak_data, 
         mapping = aes(x = XPos, xend = XPos, y = 0, yend = Height, color = SampleType),
         show.legend = FALSE
      )
   }
   
   sample_type_colors <- sapply(SAMPLE_TYPE, `[[`, "color")
   
   ggplot(plot_data, aes(x = .data[[x]], y = FeatureValue, group = SampleType)) +
      
      geom_ribbon(aes(ymin = PctMin, ymax = PctMax, fill = SampleType), alpha = 0.1) +
      geom_line(aes(color = SampleType)) +
      geom_sample_peak +
      
      scale_color_manual(values = sample_type_colors) +
      scale_fill_manual(values = sample_type_colors) +
      
      gg_hline + 
      gg_scale_y_continuous +
      plot_labels +
      theme(legend.position = "none")
}


PLOTS[[FEATURE_TYPE$COVERAGE_DISTRIBUTION]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$COVERAGE_DISTRIBUTION)
   
   plot_labels <- labs(title = "Coverage", x = "Coverage", y = "Prop. of bases")
   
   plot_distribution(
      plot_data, x = "ReadDepth", plot_labels = plot_labels, 
      mark_sample_peak = TRUE, invert_normal = TRUE
   )
})

PLOTS[[FEATURE_TYPE$FRAG_LENGTH_DISTRIBUTION]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$FRAG_LENGTH_DISTRIBUTION)
   
   plot_labels <- labs(title = "Fragment length", x = "Fragment length", y = "Prop. of fragments")
   
   plot_distribution(
      plot_data, x = "FragLength", plot_labels = plot_labels, 
      mark_sample_peak = TRUE, invert_normal = TRUE
   )
})

PLOTS[[FEATURE_TYPE$GC_BIAS]] <- local({
   
   plot_data <- get_plot_data(FEATURE_TYPE$GC_BIAS)
   
   plot_labels <- labs(title = "GC bias", x = "GC percentage", y = "Read depth")
   
   plot_distribution(
      plot_data, x = "GCBucket", plot_labels = plot_labels, 
      mark_sample_peak = FALSE, invert_normal = FALSE
   )
})

## =============================
## Summary table
## =============================

SUMMARY_TABLE_DATA <- get_plot_data(FEATURE_TYPE$SUMMARY_TABLE)

plot_sub_table <- function(feature_group, number_format = "NUMBER", show_title = TRUE, show_sample_type_label = TRUE){

   if(FALSE){
      show_title = TRUE
      show_sample_type_label = TRUE
      feature_group = "Mutational burden"; number_format = "LOG"
      feature_group = "Contamination"; number_format = "PERCENT"
   }

   ## Prep data =============================
   plot_data <- SUMMARY_TABLE_DATA %>%
      dplyr::filter(FeatureGroup == feature_group & NumberFormat == number_format) %>%
      dplyr::mutate(
         PlotLabel = factor(PlotLabel, rev(levels(PlotLabel))),
         SampleType = factor(SampleType, rev(levels(SampleType)))
      )
   
   n_rows <- plot_data$PlotLabel %>% unique() %>% length()
   
   if(n_rows == 0)
      return(NULL)

   ## Plot config =============================
   axis_trans <- "identity"
   axis_breaks <- waiver()
   
   axis_limits <- c(
      min(c(0, plot_data$PctLimLower)), 
      max(plot_data$PctLimUpper)
   )
   
   axis_fmt_func <- function(x) return(x)
   value_fmt_func <- function(x) ifelse(x < 1, signif(x, 2), round(x, 1))

   if(number_format == "PERCENT"){
      value_fmt_func <- scales::label_percent(accuracy = 0.1)
      axis_fmt_func <- scales::label_percent()
      axis_limits <- c(0, 1)
   } else if(number_format == "LOG"){
      axis_trans <- "log1p"
      axis_breaks <- 10^(0:10)
   }

   gg_div_lines <- function(direction = "horizontal"){
      
      positions <- seq(1.5, n_rows)
      color <- if(n_rows > 1) "grey70" else "white"
      linetype <- "dotted"
      
      if(direction == "horizontal"){
         geom_hline(yintercept = positions, color = color, linetype = linetype)
      } else if(direction == "vertical"){
         geom_vline(xintercept = positions, color = color, linetype = linetype)
      } else {
         LOGGER$error("`direction` must be 'horizontal' or 'vertical'")
      }
   }
   
   ## Feature values =============================
   plot_data <- plot_data %>% dplyr::mutate(
      ValueLabel = value_fmt_func(FeatureValue),
      
      QcStatus = as.character(QcStatus),
      QcStatusEnum = case_when(
         startsWith(QcStatus, QC_STATUS$FAIL$name) ~ QC_STATUS$FAIL$name,
         startsWith(QcStatus, QC_STATUS$WARN$name) ~ QC_STATUS$WARN$name,
         .default = QC_STATUS$PASS$name
      ),
      QcStatusEnum = factor(QcStatusEnum, sapply(QC_STATUS, `[[`, "name"))
   )

   sample_type_colors <- sapply(SAMPLE_TYPE, `[[`, "color")
   qc_status_colors <- sapply(QC_STATUS, `[[`, "color")
   
   subplot_values <- ggplot(plot_data, aes(y = PlotLabel, x = SampleType, group = SampleType, color = SampleType)) +
      
      geom_label(
         aes(label = ValueLabel, fill = QcStatusEnum), 
         size = 3, label.padding = unit(4, "pt"),
         border.colour = ifelse(plot_data$QcStatusEnum == QC_STATUS$PASS$name, "#FFFFFF00", "black")
      ) +
      scale_color_manual(values = sample_type_colors) +
      scale_fill_manual(values = qc_status_colors) +
      
      gg_div_lines("horizontal") +
      
      scale_x_discrete(position = "bottom", limits = rev, drop = FALSE) +
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
   subplot_boxplot <- plot_boxplot(plot_data, x = "PlotLabel", y = "FeatureValue") +
      gg_div_lines("vertical") + 
      scale_y_continuous(
         trans = axis_trans,
         breaks = axis_breaks,
         limits = axis_limits,
         label = axis_fmt_func
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
      )

   ## Combine plots =============================
   subplots_combined <- patchwork::wrap_plots(subplot_values, subplot_boxplot, nrow = 1)
   
   subplots_combined$height <- n_rows
   
   subplots_combined
}

PLOTS[[FEATURE_TYPE$SUMMARY_TABLE]] <- local({
   
   plots <- list()
   
   plot_index <- 0
   
   for(feature_group in SUMMARY_TABLE_DATA$FeatureGroup %>% unique()){
      
      is_first_group_plot <- TRUE
      
      for(number_format in NUMBER_FORMATS){
         subplot <- plot_sub_table(feature_group = feature_group, number_format = number_format)
         
         if(is.null(subplot))
            next
         
         LOGGER$debug(" featureGroup(%s) numberFormat(%s)", feature_group, number_format)
         
         plot_index <- plot_index + 1
         
         plots[[plot_index]] <- plot_sub_table(
            feature_group = feature_group, 
            number_format = number_format,
            show_title = is_first_group_plot,
            show_sample_type_label = FALSE
         )
         
         is_first_group_plot <- FALSE
      }
   }
   
   ## Show sample type only in the last plot
   plots[[plot_index]] <- plot_sub_table(
      feature_group = feature_group, 
      number_format = number_format,
      show_title = is_first_group_plot,
      show_sample_type_label = TRUE
   )
   
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
         paste(SAMPLE_TYPE$TUMOR$human_readable_name, "sample"),
         paste(SAMPLE_TYPE$TUMOR$human_readable_name, "cohort"),
         paste(SAMPLE_TYPE$NORMAL$human_readable_name, "sample"),
         paste(SAMPLE_TYPE$NORMAL$human_readable_name, "cohort")
      ),
      
      Color = c(
         SAMPLE_TYPE$TUMOR$color,
         SAMPLE_TYPE$TUMOR$color %>% adjustcolor(alpha.f = 0.2),
         SAMPLE_TYPE$NORMAL$color,
         SAMPLE_TYPE$NORMAL$color %>% adjustcolor(alpha.f = 0.2)
      )
   ) %>% 
      mutate(Index = 1:dplyr::n())
   
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
   
   get_qc_string <- function(sample_type){
      
      df <- SUMMARY_TABLE_DATA %>% filter(
         SampleType == sample_type$name & 
         nchar(as.character(QcStatus)) != 0
      )
      
      if(nrow(df) == 0)
         return(character())
      
      qc_string <- sprintf("[%s: %s]", df$PlotLabel, df$QcStatus) %>% paste(collapse = ", ")
      qc_string <- paste0(sample_type$human_readable_name, " QC status: ", qc_string)
      return(qc_string)
   }
   
   qc_strings <- c(
      get_qc_string(SAMPLE_TYPE$TUMOR), 
      get_qc_string(SAMPLE_TYPE$NORMAL)
   )
   
   subtitle <- if(length(qc_strings) > 0) paste(qc_strings, collapse = "\n") else waiver()
   
   patchwork::plot_annotation(
      title = TUMOR_ID,
      subtitle = subtitle,
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
   )
})

## =============================
## Combine plots
## =============================

create_report <- local({
   
   LOGGER$info("Combining plots")
   
   plots <- lapply(PLOTS, function(p){ patchwork::free(p, "label") })

   plot_letter_name_map <- c(
      "A" = FEATURE_TYPE$SUMMARY_TABLE,
      
      "B" = FEATURE_TYPE$COVERAGE_DISTRIBUTION,
      "C" = FEATURE_TYPE$FRAG_LENGTH_DISTRIBUTION,
      "D" = FEATURE_TYPE$GC_BIAS,
      "E" = FEATURE_TYPE$DISCORDANT_FRAG_FREQ,
      "F" = FEATURE_TYPE$DUPLICATE_FREQ,
      "G" = FEATURE_TYPE$MISSED_VARIANT_LIKELIHOOD,
      
      "H" = FEATURE_TYPE$BQR_BY_ORIG_QUAL,
      "I" = FEATURE_TYPE$BQR_BY_SNV96_CONTEXT,
      "J" = FEATURE_TYPE$MS_INDEL_ERROR_RATES,
      "K" = FEATURE_TYPE$MS_INDEL_ERROR_BIAS,
      
      "L" = PLOT_NAME_LEGEND
   )

   plots <- plots[plot_letter_name_map]
   
   design <- "
      AABBCCDD
      AAEEFFGG
      AAHHIIII
      AAJJKKLL
   "
   
   plots_combined <- patchwork::wrap_plots(plots, design = design) + REPORT_TITLE
   
   LOGGER$info("Writing report to: %s", OUTPUT_PATH)
   ggsave(
      filename = OUTPUT_PATH, plot = plots_combined, 
      device = "pdf", width = 20, height = 12, units = "in"
   )
})



