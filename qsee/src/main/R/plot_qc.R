library(ggplot2)
library(dplyr)
library(tidyr)
library(ggh4x)
library(patchwork)

library(gt)
library(svglite)

################################
## Config
################################

# args <- commandArgs(TRUE)
# 
# TUMOR_ID <- args[1]
# NORMAL_ID <- args[2]
# SAMPLE_FEATURES_FILE <- args[3]
# COHORT_NAMED_PCT_FILE <- args[4]
# OUTPUT_PATH <- args[5]
# GLOBAL_LOG_LEVEL <- args[6]

TUMOR_ID <- "H00000098"
NORMAL_ID <- "H00000098-ref"
SAMPLE_FEATURES_FILE <- "/Users/lnguyen/Hartwig/experiments/wigits_qc/analysis/20250805_oa_run_hmf_samples/qsee_output/H00000098.qsee.vis.features.tsv.gz"
COHORT_NAMED_PERCENTILES_FILE <- "/Users/lnguyen/Hartwig/experiments/wigits_qc/analysis/20250805_oa_run_hmf_samples/qsee_output/COHORT.qsee.vis.named_percentiles.tsv.gz"
OUTPUT_PATH <- sprintf("/Users/lnguyen/Hartwig/experiments/wigits_qc/analysis/20250805_oa_run_hmf_samples/qsee_output/%s.qsee.vis.report.pdf", TUMOR_ID)
GLOBAL_LOG_LEVEL <- "DEBUG"

SAMPLE_DATA <- read.delim(SAMPLE_FEATURES_FILE)
COHORT_DATA <- read.delim(COHORT_NAMED_PERCENTILES_FILE)

################################
## Helper functions
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

logMessage <- function(log_level, string){
   
   current_time <- format(Sys.time(), "%H:%H:%OS3")
   
   log_message <- sprintf("%s [R] [%5s] %s", current_time, log_level$name, string)
   
   if(log_level$severity >= LOG_LEVEL[[LOG_LEVEL$ERROR$name]]$severity)
      stop(log_message)
   
   if(log_level$severity >= LOG_LEVEL[[GLOBAL_LOG_LEVEL]]$severity)
      message(log_message)
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
   SUMMARY_TABLE              = list(name = "SUMMARY_TABLE"),
   COVERAGE_DISTRIBUTION      = list(name = "COVERAGE_DISTRIBUTION"),
   FRAG_LENGTH_DISTRIBUTION   = list(name = "FRAG_LENGTH_DISTRIBUTION"),
   MISSED_VARIANT_LIKELIHOOD  = list(name = "MISSED_VARIANT_LIKELIHOOD"),
   DUPLICATE_FREQ             = list(name = "DUPLICATE_FREQ"),
   GC_BIAS                    = list(name = "GC_BIAS"),
   DISCORDANT_READ_STATS      = list(name = "DISCORDANT_READ_STATS"),
   BQR_BY_SNV96_CONTEXT       = list(name = "BQR_PER_SNV96_CONTEXT"),
   BQR_BY_ORIG_QUAL           = list(name = "BQR_PER_ORIG_QUAL"),
   MS_INDEL_ERROR_RATES       = list(name = "MS_INDEL_ERROR_RATES"),
   MS_INDEL_ERROR_BIAS        = list(name = "MS_INDEL_ERROR_BIAS")
)

################################
## Summary table
################################

draw_summary_table <- function(tumor_id = TUMOR_ID, normal_id = NORMAL_ID){
   
   if(FALSE){
      tumor_id = TUMOR_ID
      normal_id = NORMAL_ID
   }
   
   ## Get cohort data ================================
   sample_data <- SAMPLE_DATA %>% filter(FeatureType == FEATURE_TYPE$SUMMARY_TABLE$name & SampleId %in% c(tumor_id, normal_id))
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
         SampleType = to_upper_camel_case(SampleType),
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
      
      gt(
         groupname_col = "FeatureGroup", row_group_as_column = TRUE,
         rowname_col = "Metric",
      ) %>%
      
      ## Data values
      fmt_markdown(columns = sample_type_columns) %>%
      
      text_transform(
         locations = cells_body(columns = PercentileInCohort),
         fn = function(values){ return(quantile_plots) }
      ) %>%
      
      sub_missing(columns = sample_type_columns, missing_text = "-") %>%
      
      cols_align(columns = sample_type_columns, align = "center") %>%
      cols_align(columns = Metric, align = "right") %>%
      
      ## Header
      tab_header(
         title = tumor_id, 
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
            cells_column_labels()
         )
      ) %>%
      
      ## Cell borders
      tab_options(
         table_body.hlines.style = "none",
         table.border.top.style = "none",
         column_labels.border.top.style = "none",
         table.border.bottom.style = "none",
         table_body.border.bottom.style = "none"
         
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


################################
## Plot functions
################################

get_prelim_plot_data <- function(feature_type, tumor_id = TUMOR_ID, normal_id = NORMAL_ID){
   
   if(FALSE){
      feature_type="COVERAGE_DISTRIBUTION"
   }
   
   cohort_data <- COHORT_DATA %>% filter(FeatureType == feature_type)
   sample_data <- SAMPLE_DATA %>% filter(FeatureType == feature_type & SampleId %in% c(tumor_id, normal_id))
   
   cohort_data$GroupType <- GROUP_TYPE$COHORT$name
   sample_data$GroupType <- GROUP_TYPE$SAMPLE$name
   
   cohort_data$SampleGroup <- paste0(cohort_data$SampleType, "_", cohort_data$GroupType)
   sample_data$SampleGroup <- paste0(sample_data$SampleType, "_", sample_data$GroupType)
   
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

plot_coverage_distribution <- function(){
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$COVERAGE_DISTRIBUTION$name)
   plot_data <- plot_data %>% rename(AxisX = ReadDepth)
   
   p <- plot_distribution(plot_data, show_median = TRUE, show_cohort_lines = FALSE, invert_normal = TRUE)
   p + labs(title = "Coverage", x = "Coverage", y = "Prop. of bases")
}

plot_frag_length_distribution <- function(){
   plot_data <- get_prelim_plot_data(feature_type = FEATURE_TYPE$FRAG_LENGTH_DISTRIBUTION$name)
   plot_data <- plot_data %>% rename(AxisX = FragLength)
   
   p <- plot_distribution(plot_data, show_median = TRUE, show_cohort_lines = FALSE, invert_normal = TRUE)
   p + labs(title = "Fragment length", x = "Fragment length", y = "Prop. of fragments")
}

plot_gc_bias <- function(){
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

plot_duplicate_freq <- function(){
   
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

plot_missed_variant_likelihood <- function(){
   
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

plot_bqr_by_orig_qual <- function(){
   
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

plot_bqr_by_snv96_context <- function(){
   
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

plot_ms_indel_error_rates <- function(){
   
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

plot_msi_indel_error_bias <- function(){
   
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
## Discordant read stats 
## =============================

plot_discordant_read_stats <- function(){
   
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

PLOTS <- list()

PLOTS[[FEATURE_TYPE$SUMMARY_TABLE$name]] <- draw_summary_table() %>% gt_grob() %>% wrap_elements(full = .)

PLOTS[[FEATURE_TYPE$COVERAGE_DISTRIBUTION$name]] <- plot_coverage_distribution()
PLOTS[[FEATURE_TYPE$FRAG_LENGTH_DISTRIBUTION$name]] <- plot_frag_length_distribution()
PLOTS[[FEATURE_TYPE$GC_BIAS$name]] <- plot_gc_bias()
PLOTS[[FEATURE_TYPE$DISCORDANT_READ_STATS$name]] <- plot_discordant_read_stats()

PLOTS[[FEATURE_TYPE$DUPLICATE_FREQ$name]] <- plot_duplicate_freq()
PLOTS[[FEATURE_TYPE$MISSED_VARIANT_LIKELIHOOD$name]] <- plot_missed_variant_likelihood()

PLOTS[[FEATURE_TYPE$BQR_BY_ORIG_QUAL$name]] <- plot_bqr_by_orig_qual()
PLOTS[[FEATURE_TYPE$BQR_BY_SNV96_CONTEXT$name]] <- plot_bqr_by_snv96_context()

PLOTS[[FEATURE_TYPE$MS_INDEL_ERROR_RATES$name]] <- plot_ms_indel_error_rates()
PLOTS[[FEATURE_TYPE$MS_INDEL_ERROR_BIAS$name]] <- plot_msi_indel_error_bias()

p_combined <- 
   patchwork::wrap_plots(
      PLOTS, guides="collect", ncol = 2,
      design = "
      AAAABB
      AAAACC
      AAAADD
      AAAAEE
      FFFGGG
      HHHIII
      JJJKKK
      "
   ) & 
   theme(
      plot.margin = unit(c(t=18, r=12, b=0, l=12), "pt"),
      legend.position = "bottom",
      legend.direction = "vertical"
   )

ggsave(
   filename = OUTPUT_PATH, plot = p_combined, 
   device = "pdf", width = 14, height = 17, units = "in"
)



