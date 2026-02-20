suppressPackageStartupMessages(library(dplyr))

## Pre-processing
library(dplyr)
library(tidyr)

## Plotting
library(ggplot2)
theme_set(
   theme_bw() +
   theme(
      panel.grid = element_blank()
   )
)

library(patchwork)
library(scales)
library(ggh4x)

################################
## Config
################################

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

################################
## Logging
################################

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

################################
## Helper functions
################################

preordered_factor <- function(x){ factor(x, unique(x)) }

preordered_factors <- function(df){ 
   as.data.frame(lapply(df, function(x){
      if(is.character(x)){
         x <- factor(x, unique(x))
      }
      return(x)
   }))
}

named_vector <- function(values, names){ setNames(values, names) }

named_vector_from_df <- function(df, column){ setNames(df[[column]], rownames(df)) }

df_to_strings <- function(df, key_value_sep = "=", group_sep = ";"){
   
   key_value_strings <- lapply(colnames(df), function(colname){
      paste0(colname, key_value_sep, df[[colname]])
   })
   
   strings <- do.call(paste, c(key_value_strings, sep=group_sep))
   
   return(strings)
}

strings_to_df <- function(strings, key_value_sep = "=", group_sep = ";"){
   
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

################################
## Constants
################################

SAMPLE_TYPE <- list(
   TUMOR = list(name = "TUMOR", human_readable_name = "Tumor", color = "#D1392C"),
   NORMAL = list(name = "NORMAL", human_readable_name = "Normal", color = "#4A7DB4")
)

GROUP_TYPE <- list(
   COHORT = list(name = "COHORT"),
   SAMPLE = list(name = "SAMPLE")
)

SAMPLE_GROUP <- list(
   TUMOR_COHORT = list(name = "TUMOR_COHORT"),
   NORMAL_COHORT = list(name = "NORMAL_COHORT"),
   TUMOR_SAMPLE = list(name = "TUMOR_SAMPLE"),
   NORMAL_SAMPLE = list(name = "NORMAL_SAMPLE")
)

FEATURE_TYPE <- list(
   ## Plot functions are defined later
   SUMMARY_TABLE              = list(name = "SUMMARY_TABLE", plot_func = NULL),
   COVERAGE_DISTRIBUTION      = list(name = "COVERAGE_DISTRIBUTION", plot_func = NULL),
   FRAG_LENGTH_DISTRIBUTION   = list(name = "FRAG_LENGTH_DISTRIBUTION", plot_func = NULL),
   GC_BIAS                    = list(name = "GC_BIAS", plot_func = NULL),
   DUPLICATE_FREQ             = list(name = "DUPLICATE_FREQ", plot_func = NULL),
   DISCORDANT_FRAG_FREQ       = list(name = "DISCORDANT_FRAG_FREQ", plot_func = NULL),
   MISSED_VARIANT_LIKELIHOOD  = list(name = "MISSED_VARIANT_LIKELIHOOD", plot_func = NULL),
   BQR_BY_SNV96_CONTEXT       = list(name = "BQR_PER_SNV96_CONTEXT", plot_func = NULL),
   BQR_BY_ORIG_QUAL           = list(name = "BQR_PER_ORIG_QUAL", plot_func = NULL),
   MS_INDEL_ERROR_RATES       = list(name = "MS_INDEL_ERROR_RATES", plot_func = NULL),
   MS_INDEL_ERROR_BIAS        = list(name = "MS_INDEL_ERROR_BIAS", plot_func = NULL)
)

PERCENTILE_PREFIX <- "Pct_"

NAMED_PERCENTILES <- list(
   MIN   = list(name = "PctMin"  , colname = paste0(PERCENTILE_PREFIX,  5)),
   LOWER = list(name = "PctLower", colname = paste0(PERCENTILE_PREFIX, 25)),
   MID   = list(name = "PctMid"  , colname = paste0(PERCENTILE_PREFIX, 50)),
   UPPER = list(name = "PctUpper", colname = paste0(PERCENTILE_PREFIX, 75)),
   MAX   = list(name = "PctMax"  , colname = paste0(PERCENTILE_PREFIX, 95))
)

NUMBER_FORMATS <- list(
   NUMBER = "NUMBER",
   PERCENT = "PERCENT",
   LOG = "LOG"
)

################################
## Load data
################################

load_cohort_percentiles <- function(){
   
   LOGGER$info("Loading cohort percentiles from: %s", COHORT_PERCENTILES_FILE)
   
   cohort_percentiles <- read.delim(COHORT_PERCENTILES_FILE, na.strings = c("NA", "null"))
   
   cohort_named_percentiles <- cohort_percentiles[ sapply(NAMED_PERCENTILES, `[[`, "colname") ]
   colnames(cohort_named_percentiles) <- sapply(NAMED_PERCENTILES, `[[`, "name")
   
   cohort_named_percentiles <- cbind(
      cohort_percentiles %>% select(!starts_with(PERCENTILE_PREFIX)),
      cohort_named_percentiles
   )
   
   return(cohort_named_percentiles)
}

load_sample_features <- function(){
   LOGGER$info("Loading sample features from: %s", SAMPLE_FEATURES_FILE)
   
   sample_features <- read.delim(SAMPLE_FEATURES_FILE, na.strings = c("NA", "null"))
   sample_features <- sample_features %>% filter(SampleId %in% c(TUMOR_ID, NORMAL_ID))
   
   return(sample_features)
}

COHORT_DATA <- load_cohort_percentiles()
SAMPLE_DATA <- load_sample_features()

################################
## Plot functions
################################

get_prelim_plot_data <- function(feature_type, merge_strategy = "rbind"){
   
   ## Select rows
   cohort_data <- COHORT_DATA %>% filter(FeatureType == feature_type)
   sample_data <- SAMPLE_DATA %>% filter(FeatureType == feature_type)
   
   if(nrow(sample_data) == 0){
      return(data.frame())
   }
   
   ## Only select the features cohort data that are present in the sample
   cohort_data <- cohort_data[
      paste(cohort_data$SampleType, cohort_data$FeatureName) %in% 
      paste(sample_data$SampleType, sample_data$FeatureName)
   ,]
   
   ## Assign groupings
   cohort_data$GroupType <- GROUP_TYPE$COHORT$name
   cohort_data$SampleGroup <- paste0(cohort_data$SampleType, "_", cohort_data$GroupType)
   
   sample_data$GroupType <- GROUP_TYPE$SAMPLE$name
   sample_data$SampleGroup <- paste0(sample_data$SampleType, "_", sample_data$GroupType)
   
   
   sample_data <- sample_data %>% rename(PctMid = FeatureValue)
   
   merged_data <- bind_rows(cohort_data, sample_data)
   merged_data <- merged_data %>% select(SampleGroup, GroupType, SampleType, everything())
   
   ## Split feature names into columns
   has_multiplex_feature_names <- grepl(".+=.+", merged_data$FeatureName[1])
   if(has_multiplex_feature_names){
      merged_data <- data.frame(merged_data, strings_to_df(merged_data$FeatureName))
      merged_data$FeatureName <- NULL
   }
   
   return(merged_data)
}

plot_missing_data <- function(plot_labels = labs()){
   ggplot() +
      annotate("text", x = 0, y = 0, label = "Missing sample data") +
      plot_labels +
      theme(
         axis.text = element_blank(),
         axis.ticks = element_blank()
      )
}

## =============================
## Summary table
## =============================

get_summary_table_data <- function(){
   
   cohort_data <- COHORT_DATA %>% filter(FeatureType == FEATURE_TYPE$SUMMARY_TABLE$name)
   sample_data <- SAMPLE_DATA %>% filter(FeatureType == FEATURE_TYPE$SUMMARY_TABLE$name)
   
   merged_data <- merge(
      sample_data, cohort_data, 
      by=c("SampleType", "SourceTool", "FeatureType", "FeatureName"), 
      all = TRUE, sort = FALSE
   )
   
   merged_data <- cbind(
      merged_data,
      strings_to_df(merged_data$PlotMetadata)
   )
   merged_data$PlotMetadata <- NULL
   
   merged_data <- preordered_factors(merged_data)
   
   return(merged_data)
}

SUMMARY_TABLE_DATA <- get_summary_table_data()

plot_sub_table <- function(feature_group, number_format = "NUMBER", show_title = TRUE, show_sample_type_label = TRUE){

   if(FALSE){
      show_title = TRUE
      show_sample_type_label = TRUE

      feature_group = "Mutational burden"; number_format = "LOG"
   }

   ## Prep data =============================
   plot_data <- SUMMARY_TABLE_DATA %>%
      filter(FeatureGroup == feature_group & NumberFormat == number_format) %>%
      mutate(
         PlotLabel = factor(PlotLabel, rev(levels(PlotLabel))),
         SampleType = factor(SampleType, rev(levels(SampleType)))
      )
   
   n_rows <- plot_data$PlotLabel %>% unique() %>% length()
   
   if(n_rows == 0)
      return(NULL)

   ## Plot config =============================

   axis_trans <- "identity"
   axis_breaks <- waiver()
   axis_limits <- c(0, NA)
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

   ## Constants used multiple times =============================
   H_DIV_LINES <- geom_hline(
      yintercept = seq(1.5, n_rows),
      color = if(n_rows > 1) "grey90" else "white",
      linetype = "dotted"
   )

   X_AXIS_POSITION <- "bottom"
   POSITION_DODGE_WIDTH <- 0.5
   BOXPLOT_WIDTH <- 0.25 * (plot_data$SampleType %>% unique() %>% length())

   SAMPLE_TYPE_COLORS <- c(TUMOR = SAMPLE_TYPE$TUMOR$color, NORMAL = SAMPLE_TYPE$NORMAL$color)

   ## Feature values =============================
   plot_data <- plot_data %>%
      mutate(
         HasQcStatus = nchar(as.character(QcStatus)) > 0,
         ValueLabel = value_fmt_func(FeatureValue),
         ValueLabel = ifelse(
            HasQcStatus, paste0(ValueLabel,"\n", QcStatus),
            ValueLabel
         )
      )

   subplot_values <- ggplot(plot_data, aes(y = PlotLabel, x = SampleType, group = SampleType, color = SampleType)) +

      H_DIV_LINES +
      geom_text(aes(label = ValueLabel), size = 3, lineheight = 0.8, show.legend = FALSE) +
      scale_color_manual(values = SAMPLE_TYPE_COLORS, drop = FALSE) +
      scale_x_discrete(position = X_AXIS_POSITION, limits = rev, drop = FALSE) +
      labs(title = feature_group) +
      
      theme(
         plot.title = if(show_title) element_text() else element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.ticks.x = element_blank(),
         axis.text.x = if(show_sample_type_label) element_text() else element_blank(),
         plot.margin = margin(r = 0, b = 10)
      )

   ## Sample vs cohort =============================
   subplot_boxplot <- ggplot(plot_data, aes(y = PlotLabel, x = FeatureValue, fill = SampleType)) +

      H_DIV_LINES +
      geom_boxplot(
         aes(xmin = PctMin, xlower = PctLower, xmiddle = PctMid, xupper = PctUpper, xmax = PctMax),
         position = position_dodge(width = POSITION_DODGE_WIDTH),
         stat = "identity", width = BOXPLOT_WIDTH, alpha = 0.3, size = 0.25, color = "grey70",
         show.legend = FALSE
      ) +
      geom_point(
         position = position_dodge(width = POSITION_DODGE_WIDTH),
         shape = 21, size = 2, show.legend = FALSE
      ) +

      scale_fill_manual(values = SAMPLE_TYPE_COLORS, drop = FALSE) +
      scale_x_continuous(
         position = X_AXIS_POSITION,
         trans = axis_trans,
         breaks = axis_breaks,
         limits = axis_limits,
         label = axis_fmt_func
      ) +

      theme(
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.ticks.y = element_blank(),
         axis.text.y = element_blank(),
         axis.line.y = element_blank(),
         plot.margin = margin(r = 10, b = 10, l = 0)
      )

   ## Combine plots =============================
   subplots_combined <- patchwork::wrap_plots(subplot_values, subplot_boxplot, nrow = 1)
   
   subplots_combined$height <- n_rows
   
   subplots_combined
}

FEATURE_TYPE$SUMMARY_TABLE$plot_func <- function(feature_group){
   
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
}

## =============================
## Line / PDF
## =============================

plot_distribution <- function(plot_data, invert_normal = FALSE, show_median = FALSE, show_cohort_lines = TRUE){
   
   ## Parse input ================================
   if(invert_normal){
      
      is_normal <- plot_data$SampleType == SAMPLE_TYPE$NORMAL$name
      
      plot_data <- within(plot_data, {
         PctMin   <- ifelse(is_normal, -PctMin, PctMin)
         PctLower <- ifelse(is_normal, -PctLower, PctLower)
         PctMid   <- ifelse(is_normal, -PctMid, PctMid)
         PctUpper <- ifelse(is_normal, -PctUpper, PctUpper)
         PctMax   <- ifelse(is_normal, -PctMax, PctMax)
      })
   } 
   
   ## Format categories
   plot_data$AxisX <- as.numeric(plot_data$AxisX)
   plot_data <- preordered_factors(plot_data)
   
   ## Plot ================================
   NO_COLOR <- "#FFFFFF00"
   
   PLOT_AESTHETICS <- data.frame(
      
      row.names = c(SAMPLE_GROUP$TUMOR_COHORT$name, SAMPLE_GROUP$TUMOR_SAMPLE$name, 
                    SAMPLE_GROUP$NORMAL_COHORT$name, SAMPLE_GROUP$NORMAL_SAMPLE$name),
      
      color     = c(SAMPLE_TYPE$TUMOR$color, SAMPLE_TYPE$TUMOR$color, 
                    SAMPLE_TYPE$NORMAL$color, SAMPLE_TYPE$NORMAL$color),
      
      fill      = c(SAMPLE_TYPE$TUMOR$color, NO_COLOR, 
                    SAMPLE_TYPE$NORMAL$color, NO_COLOR),
      
      linetype  = c("11", "solid",
                    "11", "solid")
   )
   
   if(!show_cohort_lines){
      PLOT_AESTHETICS[
         c(SAMPLE_GROUP$TUMOR_COHORT$name, SAMPLE_GROUP$NORMAL_COHORT$name), 
         "color"
      ] <- NO_COLOR
   }
   
   p <- ggplot(plot_data, aes(x=AxisX, group=SampleGroup)) +
      
      geom_ribbon(aes(ymin = PctMin, ymax = PctMax, fill = SampleGroup), alpha=0.1) +
      geom_line(aes(y = PctMid, color = SampleGroup, linetype = SampleGroup)) +
      
      scale_color_manual(values = named_vector_from_df(PLOT_AESTHETICS, "color") ) +
      scale_fill_manual(values = named_vector_from_df(PLOT_AESTHETICS, "fill") ) +
      scale_linetype_manual(values = named_vector_from_df(PLOT_AESTHETICS, "linetype") )
   
   ## Customisations ================================
   if(show_median){
      median_data <- plot_data %>% 
         group_by(SampleGroup) %>%
         slice(floor(n() * 0.05):ceiling(n() * 0.95)) %>% ## Remove the bottom 5% since the coverage plot has a spike there
         reframe(
            AxisX = AxisX[which.max(abs(PctMid))],
            PctMid = PctMid[which.max(abs(PctMid))]
         )
      
      p <- p +
         geom_segment(
            data = median_data, 
            mapping = aes(
               x = AxisX, xend = AxisX, y = 0, yend = PctMid, 
               linetype = SampleGroup, color = SampleGroup
            ),
            show.legend = FALSE
         )
   }
   
   if(invert_normal){ 
      p <- p +
         geom_hline(yintercept = 0, linewidth = 0.25) +
         scale_y_continuous(labels = function(x) abs(x))
   }
   
   return(p)
}

FEATURE_TYPE$COVERAGE_DISTRIBUTION$plot_func <- function(){
   
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$COVERAGE_DISTRIBUTION$name)
   plot_data <- plot_data %>% rename(AxisX = ReadDepth)
   
   p <- plot_distribution(plot_data, show_median = TRUE, show_cohort_lines = FALSE, invert_normal = TRUE)
   p + labs(title = "Coverage", x = "Coverage", y = "Prop. of bases")
}

FEATURE_TYPE$FRAG_LENGTH_DISTRIBUTION$plot_func <- function(){
   
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$FRAG_LENGTH_DISTRIBUTION$name)
   plot_data <- plot_data %>% rename(AxisX = FragLength)
   
   p <- plot_distribution(plot_data, show_median = TRUE, show_cohort_lines = FALSE, invert_normal = TRUE)
   p + labs(title = "Fragment length", x = "Fragment length", y = "Prop. of fragments")
}

FEATURE_TYPE$GC_BIAS$plot_func <- function(){
   
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$GC_BIAS$name)
   plot_data <- plot_data %>% rename(AxisX = GCBucket)
   
   p <- plot_distribution(plot_data, show_cohort_lines = FALSE)
   p + labs(title = "GC bias", x = "GC percentage", y = "Read depth")
}

## =============================
## Dot plots
## =============================

plot_dotplot <- function(plot_data, point_size = 1, linerange_size = 0.3, hlines = NULL, vlines = NULL, facet_scales = "free_x"){
   
   REQUIRED_COLUMNS <- c("SampleGroup", "PctMin", "PctMax", "AxisX")
   missing_columns <- REQUIRED_COLUMNS[ !(REQUIRED_COLUMNS %in% colnames(plot_data)) ]
   if(length(missing_columns) > 0) stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
   
   PLOT_AESTHETICS <- data.frame(
      
      row.names = c(SAMPLE_GROUP$TUMOR_COHORT$name, SAMPLE_GROUP$TUMOR_SAMPLE$name, 
                    SAMPLE_GROUP$NORMAL_COHORT$name, SAMPLE_GROUP$NORMAL_SAMPLE$name),
      
      color     = c(SAMPLE_TYPE$TUMOR$color, SAMPLE_TYPE$TUMOR$color,
                    SAMPLE_TYPE$NORMAL$color, SAMPLE_TYPE$NORMAL$color),
      
      linetype  = c("11", "blank",
                    "11", "blank"),
      
      alpha     = c(0, 1,
                    0, 1)
   )
   
   DODGE_WIDTH <- 0.5
   
   should_facet_x <- "FacetX" %in% colnames(plot_data)
   should_facet_y <- "FacetY" %in% colnames(plot_data)
   facet_formula <- 
      if(should_facet_x & should_facet_y){ 
         "FacetY ~ FacetX" 
      } else if(should_facet_x){ 
         ". ~ FacetX" 
      } else if(should_facet_y){ 
         "FacetY ~ ." 
      } else { NULL }
   
   p <- ggplot(plot_data, aes(x = AxisX, y = PctMid, group = SampleType)) +
      
      { if(!is.null(facet_formula)) facet_grid(facet_formula, scales = facet_scales) } +
      
      { if(!is.null(hlines)) geom_hline(linewidth = 0.25, color = "grey", yintercept = hlines) } +
      { if(!is.null(vlines)) geom_vline(linewidth = 0.25, color = "grey", xintercept = vlines) } +
      
      geom_linerange(
         aes(ymin = PctMin, ymax = PctMax, color = SampleGroup, linetype = SampleGroup),
         position = position_dodge(width = DODGE_WIDTH), linewidth = 0.3
      ) +
      scale_linetype_manual(values = named_vector_from_df(PLOT_AESTHETICS, "linetype") ) +
      
      geom_point(
         aes(color = SampleGroup, alpha = SampleGroup),
         position = position_dodge(width = DODGE_WIDTH), size = 1
      ) +
      scale_color_manual(values = named_vector_from_df(PLOT_AESTHETICS, "color") ) +
      scale_alpha_manual(values = named_vector_from_df(PLOT_AESTHETICS, "alpha") ) +
      
      theme(
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         panel.spacing.x = unit(-0.5, "pt"),
         axis.text.y.right = element_blank(),
         axis.ticks.y.right = element_blank(),
      )
   
   return(p)
}

FEATURE_TYPE$BQR_BY_ORIG_QUAL$plot_func <- function(){
   
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$BQR_BY_ORIG_QUAL$name)
   
   plot_labels <- labs(
      title = "BQR by original base quality", 
      x = "Original base quality", 
      y = "Phred score adjustment"
   )
   
   if(nrow(plot_data) == 0){
      return(plot_missing_data(plot_labels))
   }
   
   plot_data <- plot_data %>% rename(
      AxisX = OriginalQualBin, 
      FacetX = StandardMutation, 
      FacetY = ReadType
   )
   
   plot_data <- preordered_factors(plot_data)
   
   plot_dotplot(plot_data, hlines = 0, facet_scales = "free") + 
      scale_y_continuous(
         labels = function(x){ ifelse(x > 0, paste0("+",x), x) },
         sec.axis = dup_axis(name = "Consensus type")
      ) +
      
      plot_labels +
      theme(
         axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
         axis.text.y.right = element_blank(),
         axis.ticks.y.right = element_blank(),
      )
}

FEATURE_TYPE$BQR_BY_SNV96_CONTEXT$plot_func <- function(){
   
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$BQR_BY_SNV96_CONTEXT$name)
   
   plot_labels <- labs(
      title = if(nrow(plot_data) == 0){
         "BQR by SNV96 context"
      } else {
         paste0("BQR by SNV96 context, base quality: ", plot_data$OriginalQualBin[1])
      },
      x = "Mutation context", 
      y = "Phred score adjustment"
   )
   
   if(nrow(plot_data) == 0){
      return(plot_missing_data(plot_labels))
   }
   
   plot_data <- plot_data %>% rename(
      AxisX = StandardTrinucContext, 
      FacetX = StandardMutation, 
      FacetY = ReadType
   )
   
   plot_data <- preordered_factors(plot_data)
   
   plot_dotplot(plot_data, hlines = 0, vlines = c(4, 8, 12) + 0.5) + 
      
      scale_y_continuous(
         labels = function(x){ ifelse(x > 0, paste0("+",x), x) },
         sec.axis = dup_axis(name = "Consensus type")
      ) +
      
      plot_labels +
      theme(
         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
         axis.text.y.right = element_blank(),
         axis.ticks.y.right = element_blank()
      )
}

FEATURE_TYPE$MS_INDEL_ERROR_RATES$plot_func <- function(){
   
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$MS_INDEL_ERROR_RATES$name)
   
   plot_labels <- labs(title = "Microsatellite indel error rates", x = "Repeat units", y = "Phred score") 
   
   if(nrow(plot_data) == 0){
      return(plot_missing_data(plot_labels))
   }
   
   plot_data <- plot_data %>% rename(
      AxisX = RefNumUnits, 
      FacetX = RepeatUnitType, 
      FacetY = ConsensusType
   )
   
   plot_data <- preordered_factors(plot_data)
   
   plot_dotplot(plot_data) +
      scale_x_discrete(breaks = function(x) ifelse(as.numeric(x) %% 3 == 0, x, "") ) +
      scale_y_continuous(limits = c(0, NA), sec.axis = dup_axis(name = "Consensus type")) +
      plot_labels +
      theme(
         panel.grid.major = element_line(color = "grey90", linewidth = 0.25)
      )
}

FEATURE_TYPE$MS_INDEL_ERROR_BIAS$plot_func <- function(){
   
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$MS_INDEL_ERROR_BIAS$name)
   
   plot_labels <- labs(title = "Microsatellite indel error bias", x = "Repeat units", y = "Phred score diff.")
   
   if(nrow(plot_data) == 0){
      return(plot_missing_data(plot_labels))
   }
   
   plot_data <- plot_data %>% rename(
      AxisX = RefNumUnits, 
      FacetX = RepeatUnitType, 
      FacetY = ConsensusType
   )
   
   plot_data <- preordered_factors(plot_data)
   
   text_data <- data.frame(
      SampleType = SAMPLE_TYPE$TUMOR$name,
      FacetX     = plot_data$FacetX[1],
      AxisX      = c(-Inf, -Inf),
      PctMid     = c( Inf, -Inf),
      Label      = c("More ins. errors", "More del. errors"),
      Vjust      = c(1.5, -0.5)
   )
   
   plot_dotplot(plot_data, data_type, hlines = 0) +
      
      geom_text(
         data = text_data,
         mapping = aes(x = AxisX, y = PctMid, label = Label, vjust = Vjust),
         hjust = -0.05, size = 2.5
      ) +
      
      scale_x_discrete(breaks = function(x) ifelse(as.numeric(x) %% 3 == 0, x, "") ) +
      
      scale_y_continuous(
         labels = function(x){ ifelse(x > 0, paste0("+",x), x) },
         sec.axis = dup_axis(name = "Consensus type")
      ) +
      
      plot_labels
}

## =============================
## Box plots
## =============================

plot_boxplot <- function(plot_data){
   
   REQUIRED_COLUMNS <- c("SampleGroup", "PctMin", "PctLower", "PctMid", "PctUpper", "PctMax", "AxisX")
   missing_columns <- REQUIRED_COLUMNS[ !(REQUIRED_COLUMNS %in% colnames(plot_data)) ]
   if(length(missing_columns) > 0) stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
   
   PLOT_AESTHETICS <- data.frame(
      
      row.names = c(SAMPLE_GROUP$TUMOR_COHORT$name, SAMPLE_GROUP$TUMOR_SAMPLE$name, 
                    SAMPLE_GROUP$NORMAL_COHORT$name, SAMPLE_GROUP$NORMAL_SAMPLE$name),
      
      color     = c(SAMPLE_TYPE$TUMOR$color, SAMPLE_TYPE$TUMOR$color,
                    SAMPLE_TYPE$NORMAL$color, SAMPLE_TYPE$NORMAL$color)
   )
   
   ggplot(plot_data, aes(x = AxisX)) +
      
      geom_boxplot(
         data = subset(plot_data, GroupType == GROUP_TYPE$COHORT$name),
         mapping = aes(
            ymin = PctMin, lower = PctLower, middle = PctMid, upper = PctUpper, ymax = PctMax,
            fill = SampleGroup
         ),
         position = position_dodge(width = 0.5),
         stat = "identity", width = 0.5, alpha = 0.3, size = 0.25, color = "grey"
      ) +
      
      geom_point(
         data = subset(plot_data, GroupType == GROUP_TYPE$SAMPLE$name),
         mapping = aes(y = PctMid, fill = SampleGroup),
         position = position_dodge(width = 0.5),
         shape = 21
      ) +
      
      scale_color_manual(values = named_vector_from_df(PLOT_AESTHETICS, "color")) +
      scale_fill_manual(values = named_vector_from_df(PLOT_AESTHETICS, "color")) +
      
      theme(
         panel.grid.minor = element_blank(),
         panel.grid.major.x = element_blank(),
         strip.text.x = element_blank(),
         axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
         legend.position = "none"
      )
}

FEATURE_TYPE$DUPLICATE_FREQ$plot_func <- function(){
   
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$DUPLICATE_FREQ$name) %>% preordered_factors()
   
   plot_labels <- labs(title = "Duplicate frequency", x = "Duplicate read count", y = "Prop. of read groups")
   
   if(nrow(plot_data) == 0){
      return(plot_missing_data(plot_labels))
   }
   
   plot_data <- plot_data %>% rename(AxisX = ReadCount)
   
   plot_boxplot(plot_data) +
      plot_labels +
      theme(
         panel.grid.major.y = element_line(color = "grey", linewidth = 0.25),
         axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
      )
}

FEATURE_TYPE$MISSED_VARIANT_LIKELIHOOD$plot_func <- function(){
   
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$MISSED_VARIANT_LIKELIHOOD$name)
   
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
      filter(GroupType == GROUP_TYPE$SAMPLE$name & PctMid >= MIN_MISSED_VARIANT_LIKELIHOOD) %>%
      arrange(-PctMid) %>%
      pull(Gene) %>% 
      unique() %>% 
      head(TOP_N_GENES)
   
   plot_data <- plot_data %>% 
      filter(Gene %in% sample_genes_of_interest) %>% 
      mutate(Gene = factor(Gene, sample_genes_of_interest)) %>%
      rename(AxisX = Gene)
   
   plot_data <- preordered_factors(plot_data)
   
   plot_boxplot(plot_data) + 
      plot_labels +
      theme(
         axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      )
}

FEATURE_TYPE$DISCORDANT_FRAG_FREQ$plot_func <- function(){
   
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$DISCORDANT_FRAG_FREQ$name) %>% preordered_factors()
   
   plot_labels <- labs(title = "Discordant fragment frequency", x = "Discordant fragment type", y = "Prop. of reads")
   
   if(nrow(plot_data) == 0){
      return(plot_missing_data(plot_labels))
   }
   
   plot_data <- plot_data %>% rename(AxisX = DiscordantFragType)
   
   plot_boxplot(plot_data) + 
      plot_labels + 
      scale_y_continuous(
         transform = "log10", 
         labels = scales::trans_format("log10", scales::math_format(10^.x))
      )
}

################################
## Combine plots
################################

create_report <- function(){

   plots <- list()
   for(i in 1:length(FEATURE_TYPE)){
      feature_type <- FEATURE_TYPE[[i]]
      LOGGER$info("Plotting featureType(%s)", feature_type$name)
      plots[[feature_type$name]] <- feature_type$plot_func() %>% patchwork::free("label")
   }
   
   plots <- plots[c(
      FEATURE_TYPE$SUMMARY_TABLE$name,
      
      FEATURE_TYPE$COVERAGE_DISTRIBUTION$name,
      FEATURE_TYPE$FRAG_LENGTH_DISTRIBUTION$name,
      FEATURE_TYPE$GC_BIAS$name,
      
      FEATURE_TYPE$DISCORDANT_FRAG_FREQ$name,
      FEATURE_TYPE$DUPLICATE_FREQ$name,
      FEATURE_TYPE$MISSED_VARIANT_LIKELIHOOD$name,
      
      FEATURE_TYPE$BQR_BY_ORIG_QUAL$name,
      FEATURE_TYPE$BQR_BY_SNV96_CONTEXT$name,
      FEATURE_TYPE$MS_INDEL_ERROR_RATES$name,
      FEATURE_TYPE$MS_INDEL_ERROR_BIAS$name
   )]
   
   design <- "
      AABBEE
      AACCFF
      AADDGG
      AAHHHH
      AAIIII
   "
   
   plots_combined <- 
      patchwork::wrap_plots(plots, guides = "collect", design = design) &
      theme(
         legend.position = "bottom",
         legend.byrow = TRUE,
         legend.direction = "vertical"
      )
   
   LOGGER$info("Writing report to: %s", OUTPUT_PATH)
   ggsave(
      filename = OUTPUT_PATH, plot = plots_combined, 
      device = "pdf", width = 16, height = 13, units = "in"
   )
}

create_report()



