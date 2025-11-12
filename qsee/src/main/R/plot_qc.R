suppressPackageStartupMessages(library(dplyr))

## Pre-processing
library(dplyr)
library(tidyr)

## Plotting
library(ggplot2)
library(ggh4x)
library(patchwork)

## html/svg/png conversion
library(svglite)
library(gt)

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
    TUMOR_ID <- "H00000098"
    NORMAL_ID <- "H00000098-ref"
    SAMPLE_FEATURES_FILE <- "/Users/lnguyen/Hartwig/experiments/wigits_qc/analysis/20250805_oa_run_hmf_samples/qsee_output/H00000098.qsee.vis.features.tsv.gz"
    COHORT_PERCENTILES_FILE <- "/Users/lnguyen/Hartwig/experiments/wigits_qc/analysis/20250805_oa_run_hmf_samples/qsee_output/COHORT.qsee.percentiles.tsv.gz"
    OUTPUT_PATH <- sprintf("/Users/lnguyen/Hartwig/experiments/wigits_qc/analysis/20250805_oa_run_hmf_samples/qsee_output/%s.qsee.vis.report.pdf", TUMOR_ID)
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
LOGGER$debug(" cohort_percentiles_file:%s", COHORT_PERCENTILES_FILE)
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
   
   group_strings <- strsplit(as.character(strings), group_sep, fixed=TRUE) %>% do.call(rbind, .)
   
   column_names <- sub(paste0(key_value_sep, ".+"), "", group_strings[1,])
   value_strings <- sub(paste0(".+", key_value_sep), "", group_strings)
   
   df <- data.frame(value_strings)
   colnames(df) <- column_names
   
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
   ## Notes:
   ## - Plot functions are defined later
   ## - The order of feature types defined here determines the plot order
   SUMMARY_TABLE              = list(name = "SUMMARY_TABLE", plot_func = NULL),
   COVERAGE_DISTRIBUTION      = list(name = "COVERAGE_DISTRIBUTION", plot_func = NULL),
   FRAG_LENGTH_DISTRIBUTION   = list(name = "FRAG_LENGTH_DISTRIBUTION", plot_func = NULL),
   MISSED_VARIANT_LIKELIHOOD  = list(name = "MISSED_VARIANT_LIKELIHOOD", plot_func = NULL),
   DUPLICATE_FREQ             = list(name = "DUPLICATE_FREQ", plot_func = NULL),
   GC_BIAS                    = list(name = "GC_BIAS", plot_func = NULL),
   DISCORDANT_READ_STATS      = list(name = "DISCORDANT_READ_STATS", plot_func = NULL),
   BQR_BY_SNV96_CONTEXT       = list(name = "BQR_PER_SNV96_CONTEXT", plot_func = NULL),
   BQR_BY_ORIG_QUAL           = list(name = "BQR_PER_ORIG_QUAL", plot_func = NULL),
   MS_INDEL_ERROR_RATES       = list(name = "MS_INDEL_ERROR_RATES", plot_func = NULL),
   MS_INDEL_ERROR_BIAS        = list(name = "MS_INDEL_ERROR_BIAS", plot_func = NULL)
)

NAMED_PERCENTILES <- list(
   MIN   = list(name = "Min"  , percentile =  5.0),
   LOWER = list(name = "Lower", percentile = 25.0),
   MID   = list(name = "Mid"  , percentile = 50.0),
   UPPER = list(name = "Upper", percentile = 75.0),
   MAX   = list(name = "Max"  , percentile = 95.0)
)

################################
## Load data
################################

load_cohort_percentiles <- function(){
   
   LOGGER$info("Loading cohort percentiles from: %s", COHORT_PERCENTILES_FILE)
   
   PERCENTILE_PREFIX <- "Pct"
   
   cohort_percentiles <- read.delim(COHORT_PERCENTILES_FILE)
   
   pct_ref_values <- cohort_percentiles %>% select(starts_with(PERCENTILE_PREFIX))
   colnames(pct_ref_values) <- sub(paste0(PERCENTILE_PREFIX, "_"), "", colnames(pct_ref_values))
   percentiles <- as.numeric(colnames(pct_ref_values))

   get_pct_ref_values <- function(target_percentile){

      if(target_percentile %in% percentiles){
         output <- pct_ref_values[[as.character(target_percentile)]]
      } else { 
         ## Linear interpolation
         
         LOGGER$debug("Interpolating percentile(%s) as it was not found in the cohort percentiles file", target_percentile)
         
         lower_index <- abs(target_percentile - percentiles) %>% which.min()
         upper_index <- lower_index + 1
         
         lower_percentile <- percentiles[lower_index]
         upper_percentile <- percentiles[upper_index]
         
         lower_ref_values <- pct_ref_values[[as.character(lower_percentile)]]
         upper_ref_values <- pct_ref_values[[as.character(upper_percentile)]]
         
         fraction <- target_percentile - lower_percentile
         output <- lower_ref_values + fraction * (upper_ref_values - lower_ref_values)
      }
      
      return(output)
   }
   
   cohort_named_percentiles <- cohort_percentiles %>% select(!starts_with(PERCENTILE_PREFIX))
   
   for(named_percentile in NAMED_PERCENTILES){
      percentile_colname <- paste0(PERCENTILE_PREFIX, named_percentile$name)
      ref_values <- get_pct_ref_values(named_percentile$percentile)
      cohort_named_percentiles[[percentile_colname]] <- ref_values
   }
   
   return(cohort_named_percentiles)
}

load_sample_features <- function(){
   LOGGER$info("Loading sample features from: %s", SAMPLE_FEATURES_FILE)
   
   sample_features <- read.delim(SAMPLE_FEATURES_FILE)
   sample_features <- sample_features %>% filter(SampleId %in% c(TUMOR_ID, NORMAL_ID))
   
   return(sample_features)
}

COHORT_DATA <- load_cohort_percentiles()
SAMPLE_DATA <- load_sample_features()

################################
## Summary table
################################

draw_summary_table <- function(){
   
   ## Get cohort data ================================
   sample_data <- SAMPLE_DATA %>% filter(FeatureType == FEATURE_TYPE$SUMMARY_TABLE$name)
   cohort_data <- COHORT_DATA %>% filter(FeatureType == FEATURE_TYPE$SUMMARY_TABLE$name)
   
   rownames(sample_data) <- sample_data %>% select(SampleType, FeatureName) %>% df_to_strings()
   rownames(cohort_data) <- cohort_data %>% select(SampleType, FeatureName) %>% df_to_strings()
   
   table_data_long <- cbind(
      sample_data,
      sample_data$FeatureName %>% strings_to_df(),
      cohort_data[rownames(sample_data),] %>% select(PctMin, PctMax)
   )
   rownames(table_data_long) <- NULL
   
   ## Formatting ================================
   format_numbers <- function(numbers){
      ifelse(
         numbers == floor(numbers), 
         formatC(numbers, format="f", digits=0),
         formatC(numbers, format = "fg", digits = 2)
      ) %>% trimws()
   }
   
   SAMPLE_TYPE_COLORS <- c(TUMOR = SAMPLE_TYPE$TUMOR$color, NORMAL = SAMPLE_TYPE$NORMAL$color)
   
   table_data_long <- table_data_long %>% 
      
      mutate(
         FeatureValueMarkdown = sprintf(
            "<span style='color:%s;'>%s</span>", 
            SAMPLE_TYPE_COLORS[SampleType], 
            format_numbers(FeatureValue)
         )
      ) %>%
      
      mutate(
         FeatureValue = format_numbers(FeatureValue),
         PctMin = format_numbers(PctMin),
         PctMax = format_numbers(PctMax),
      ) %>%
      
      preordered_factors()
   
   ## Wide table ================================
   table_data_wide <- table_data_long %>% 
      select(FeatureGroup, Metric, SampleType, FeatureValueMarkdown) %>%
      mutate(
         SampleType = sapply(SampleType, function(x){ SAMPLE_TYPE[[x]]$human_readable_name }),
      ) %>%
      pivot_wider(
         names_from = SampleType,
         values_from = FeatureValueMarkdown,
         names_vary = "slowest"
      ) %>%
      as.data.frame()
   
   table_data_wide$PercentileInCohort <- NA ## Placeholder; plots will be inserted later into this column
   
   ## Plot sample quantile in cohort ================================
   DODGE_DISTANCE <- 0.1
   
   quantile_plots <- lapply(unique(table_data_wide$Metric), function(metric){
      #metric="DualStrandReadsRate"
      #metric="MeanCoverage"
      
      plot_data <- table_data_long %>% filter(Metric == metric)
      plot_data$YPosition <- 1:nrow(plot_data) * DODGE_DISTANCE
      
      p <- ggplot(plot_data, aes(x = PctInCohort, y = YPosition, fill = SampleType)) +
         
         geom_vline(
            xintercept = c(0, 25, 50, 75, 100),
            color = c("black", "lightgrey", "lightgrey", "lightgrey", "black"),
            linewidth = 0.4
         ) +
         
         geom_point(shape = 21, size = 2.5) +
         
         scale_x_continuous(limits = c(0, 100)) +
         scale_y_continuous(expand = c(DODGE_DISTANCE, DODGE_DISTANCE)) +
         scale_fill_manual(values = SAMPLE_TYPE_COLORS, guide="none") +
         
         theme_void() +
         theme(
            panel.grid = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "pt")
         )
      
      svg_string_builder <- svgstring(height = 0.2, width = 1)
      plot(p)
      dev.off()
      svg_string <- svg_string_builder()
      
      gt::html(svg_string)
   })
   
   sample_type_columns <- table_data_long$SampleType %>% 
      unique() %>% 
      sapply(., function(x){ SAMPLE_TYPE[[x]]$human_readable_name })
   
   table_data_wide %>%
      
      gt(groupname_col = "FeatureGroup", rowname_col = "Metric") %>%
      
      ## Data values
      fmt_markdown(columns = sample_type_columns) %>%
      
      text_transform(
         locations = cells_body(columns = PercentileInCohort),
         fn = function(values){ return(quantile_plots) }
      ) %>%
      
      sub_missing(columns = sample_type_columns, missing_text = "-") %>%
      
      cols_align(columns = sample_type_columns, align = "center") %>%
      cols_align(columns = Metric, align = "left") %>%
      
      ## Header
      tab_header(
         title = TUMOR_ID,
         #subtitle = paste0(names(sample_qc), ": ", sample_qc) %>% paste(collapse = "<br>") %>% md()
      ) %>%
      
      tab_options(
         heading.align = "left",
         heading.title.font.size = 24,
         heading.subtitle.font.size = 16,
      ) %>%
      
      tab_style(
         style = cell_text(weight = "bold"),
         locations = list(
            cells_title(groups = "title"),
            cells_column_labels(),
            cells_row_groups()
         )
      ) %>%
      
      ## Cell borders
      tab_options(
         table.border.top.style = "none",
         table.border.bottom.style = "none",
         table_body.border.bottom.style = "none",
         table_body.hlines.style = "none",
         column_labels.border.top.style = "none"
      ) %>%
      
      tab_style(
         style = cell_borders(sides = "left", color = "lightgrey", weight = px(1), style = "solid"),
         locations = cells_body(columns = contains("Value"))
      ) %>%
      
      tab_style(
         style = cell_borders(sides = "top", color = "lightgrey", weight = px(1), style = "solid"),
         locations = list(
            cells_row_groups(),
            cells_stub(rows = !duplicated(FeatureGroup)),
            cells_body(rows = !duplicated(FeatureGroup))
         )
      )
}

gt_grob <- function(gt_object, ...){
   
   out_name <- file.path(tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png"))
   gt::gtsave(gt_object, out_name, ...)
   
   in_png <- png::readPNG(out_name)
   on.exit(file.remove(out_name), add=TRUE)
   
   grid::rasterGrob(in_png)
}

FEATURE_TYPE$SUMMARY_TABLE$plot_func <- function(){
   draw_summary_table() %>% gt_grob() %>% wrap_elements(full = .)
}


################################
## Plot functions
################################

get_prelim_plot_data <- function(feature_type){
   
   ## Select rows
   cohort_data <- COHORT_DATA %>% filter(FeatureType == feature_type)
   sample_data <- SAMPLE_DATA %>% filter(FeatureType == feature_type)
   
   ## Assign groupings
   cohort_data$GroupType <- GROUP_TYPE$COHORT$name
   sample_data$GroupType <- GROUP_TYPE$SAMPLE$name
   
   cohort_data$SampleGroup <- paste0(cohort_data$SampleType, "_", cohort_data$GroupType)
   sample_data$SampleGroup <- paste0(sample_data$SampleType, "_", sample_data$GroupType)
   
   ## Merge cohort and sample data into one data frame
   sample_data <- sample_data %>% rename(PctMid = FeatureValue)
   
   sample_data <- sapply(colnames(cohort_data), function(column){
      if(!(column %in% colnames(sample_data))){
         return(NA)
      } else {
         return(sample_data[,column])
      }
   }) %>% as.data.frame()
   
   plot_data <- rbind(cohort_data, sample_data)
   plot_data <- plot_data %>% select(SampleGroup, GroupType, SampleType, everything())
   
   ## Split feature names into columns
   has_multiplex_feature_names <- grepl(".+=.+", plot_data$FeatureName[1])
   if(has_multiplex_feature_names){
      plot_data <- data.frame(plot_data, strings_to_df(plot_data$FeatureName))
      plot_data$FeatureName <- NULL
   }
   
   return(plot_data)
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
      scale_linetype_manual(values = named_vector_from_df(PLOT_AESTHETICS, "linetype") ) +
      
      theme_bw() +
      theme(panel.grid = element_blank())
   
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
   p + labs(title = "Read depth", x = "GC percentage", y = "Read depth")
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
      
      ##labs(title = data_type$name, x = data_type$axis_x_title, y = data_type$axis_y_title) +
      
      theme_bw() +
      theme(
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         panel.spacing.x = unit(-0.5, "pt"),
         axis.text.y.right = element_blank(),
         axis.ticks.y.right = element_blank(),
      )
   
   return(p)
}

FEATURE_TYPE$DUPLICATE_FREQ$plot_func <- function(){
   
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$DUPLICATE_FREQ$name)
   plot_data <- preordered_factors(plot_data)
   plot_data <- plot_data %>% rename(AxisX = ReadCount)
   
   plot_dotplot(plot_data) +
      
      labs(title = "Duplicate frequency", x = "Duplicate read count", y = "Prop. of read groups") +
      
      theme(
         panel.grid.major.y = element_line(),
         axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
      )
}

FEATURE_TYPE$MISSED_VARIANT_LIKELIHOOD$plot_func <- function(){
   
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$MISSED_VARIANT_LIKELIHOOD$name)
   
   MIN_MISSED_VARIANT_LIKELIHOOD <- 0.01
   TOP_N_GENES <- 20
   
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
   
   plot_dotplot(plot_data, data_type) + 
      labs(
         title = sprintf("Top %s genes with potential missed variants", TOP_N_GENES),
         x = "Gene", y = "Missed variant likelihood"
      ) +
      theme(
         axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      )
}

FEATURE_TYPE$BQR_BY_ORIG_QUAL$plot_func <- function(){
   
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$BQR_BY_ORIG_QUAL$name)
   
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
      
      labs(
         title = " BQR by original base quality", 
         x = "Original base quality", 
         y = "Phred score adjustment"
      ) +
      
      theme(
         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
         axis.text.y.right = element_blank(),
         axis.ticks.y.right = element_blank(),
      )
}

FEATURE_TYPE$BQR_BY_SNV96_CONTEXT$plot_func <- function(){
   
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$BQR_BY_SNV96_CONTEXT$name)

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
      
      labs(
         title = "BQR by SNV96 context (qual. 30+)", 
         x = "Mutation context", 
         y = "Phred score adjustment"
      ) +
      
      theme(
         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
         axis.text.y.right = element_blank(),
         axis.ticks.y.right = element_blank()
      )
}

FEATURE_TYPE$MS_INDEL_ERROR_RATES$plot_func <- function(){
   
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$MS_INDEL_ERROR_RATES$name)
   
   plot_data <- plot_data %>% rename(
      AxisX = RefNumUnits, 
      FacetX = RepeatUnitType, 
      FacetY = ConsensusType
   )
   
   plot_data <- preordered_factors(plot_data)
   
   plot_dotplot(plot_data) +
      
      scale_x_discrete(breaks = function(x) ifelse(as.numeric(x) %% 3 == 0, x, "") ) +
      scale_y_reverse(sec.axis = dup_axis(name = "Consensus type")) +
      
      labs(title = "Microsatellite indel error rates", x = "Repeat units", y = "Phred score") +
      
      theme(
         panel.grid.major = element_line()
      )
}

FEATURE_TYPE$MS_INDEL_ERROR_BIAS$plot_func <- function(){
   
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$MS_INDEL_ERROR_BIAS$name)
   
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
      
      labs(title = "Microsatellite indel error bias", x = "Repeat units", y = "Phred score diff.") +
      
      scale_y_continuous(
         labels = function(x){ ifelse(x > 0, paste0("+",x), x) },
         sec.axis = dup_axis(name = "Consensus type")
      )
}

## =============================
## Box plots
## =============================

FEATURE_TYPE$DISCORDANT_READ_STATS$plot_func <- function(){
   
   return(ggplot() + theme_bw())
   
   # plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$DISCORDANT_READ_STATS$name)
   # plot_data <- preordered_factors(plot_data)
   # 
   # PLOT_AESTHETICS <- data.frame(
   #    
   #    row.names = c(SAMPLE_GROUP$TUMOR_COHORT$name, SAMPLE_GROUP$TUMOR_SAMPLE$name,
   #                  SAMPLE_GROUP$NORMAL_COHORT$name, SAMPLE_GROUP$NORMAL_SAMPLE$name),
   #    
   #    color     = c("#D1392C66", SAMPLE_TYPE$TUMOR$color, 
   #                  "#4A7DB466", SAMPLE_TYPE$NORMAL$color)
   # )
   # 
   # ggplot(plot_data, aes(x = FeatureName)) +
   #    
   #    geom_boxplot(
   #       data = subset(plot_data, GroupType == GROUP_TYPE$COHORT$name),
   #       mapping = aes(
   #          ymin = PctMin, lower = PctLower, middle = PctMid, upper = PctUpper, ymax = PctMax, 
   #          color = SampleGroup
   #       ),
   #       stat = "identity"
   #    ) +
   #    
   #    geom_point(
   #       data = subset(plot_data, GroupType == GROUP_TYPES$sample),
   #       mapping = aes(y = QuantileMid, color = SampleGroup)
   #    ) +
   #    
   #    scale_color_manual(values = named_vector_from_df(PLOT_AESTHETICS, "color")) +
   #    scale_y_continuous(transform = "log10") +
   #    
   #    labs(title = data_type$name, x = data_type$axis_x_title, y = data_type$axis_y_title) +
   #    
   #    theme_bw() +
   #    theme(
   #       panel.grid.minor = element_blank(),
   #       panel.grid.major.x = element_blank(),
   #       strip.text.x = element_blank(),
   #       axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
   #    )
}


################################
## Combine plots
################################

create_report <- function(){

   LOGGER$info("Creating plots per feature type")

   plots <- list()

   for(feature_type in FEATURE_TYPE){
      LOGGER$debug("Plotting: %s", feature_type$name)
      plots[[feature_type$name]] <- feature_type$plot_func()
   }
   
   plots_combined <- 
      patchwork::wrap_plots(
         plots, guides="collect", ncol = 2,
         design = "
         AABB
         AACC
         AADD
         AAEE
         FFGG
         HHII
         JJKK
         "
      ) & 
      theme(
         plot.margin = unit(c(t=18, r=12, b=0, l=12), "pt"),
         legend.position = "bottom",
         legend.direction = "vertical"
      )
   
   LOGGER$info("Writing report to: %s", OUTPUT_PATH)
   ggsave(
      filename = OUTPUT_PATH, plot = plots_combined, 
      device = "pdf", width = 14, height = 18, units = "in"
   )
}

create_report()



