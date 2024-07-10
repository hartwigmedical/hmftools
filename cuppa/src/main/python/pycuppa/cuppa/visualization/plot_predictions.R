options(max.print=500)
options(stringsAsFactors=FALSE)

## Args ================================
args <- commandArgs(trailingOnly = TRUE)

VIS_DATA_PATH <- args[1]
PLOT_PATH <- args[2]

## Dependencies ================================
load_packages <- function(){
   required_packages <- c("ggplot2", "ggh4x", "stringr", "patchwork")
   
   missing_packages <- required_packages[!(required_packages %in% rownames(installed.packages()))]
   if(length(missing_packages) > 0){
      warning("Installing missing packages: ", paste(missing_packages, collapse=", "))
      install.packages(missing_packages, repos="https://cloud.r-project.org")
   }
   
   for(package in required_packages)
      suppressWarnings(suppressMessages(library(package, character.only=TRUE)))
}

load_packages()

## Load data ================================
VIS_DATA <- read.delim(VIS_DATA_PATH)

CANCER_TYPE_METADATA <- (function(){
   
   cancer_types <- unique(VIS_DATA$cancer_type)
   
   prefix <- str_extract(cancer_types, "^[A-Za-z/ ]+")
   suffix <- str_replace(cancer_types, "^[A-Za-z/ ]+: ", "")
   
   df <- data.frame(supertype = prefix, subtype = suffix, row.names = cancer_types)
   
   ## Convert character vectors to factors
   for(colname in colnames(df)){
      if(is.character(df[[colname]])){
         df[[colname]] <- factor(df[[colname]], unique(df[[colname]]))
      }
   }
   
   subtypes_per_supertype <- sapply(split(df$subtype, df$supertype), length)
   df$subtypes_per_supertype <- subtypes_per_supertype[df$supertype]
   
   return(df)
})()

VIS_DATA$cancer_supertype <- CANCER_TYPE_METADATA[VIS_DATA$cancer_type, "supertype"]

## Constants ================================
MAPPINGS_CLASSIFIER_NAMES <- c(
   "combined"="COMBINED",
   
   "dna_combined"="DNA COMBINED",
   "gen_pos"="GENOMIC POSITION",
   "snv96"="SNV96",
   "event"="EVENT",
   
   "rna_combined"="RNA COMBINED",
   "gene_exp"="GENE EXPRESSION",
   "alt_sj"="ALT SJ"
)

MAPPING_EVENT_NAMES <- c(
   "driver"="driver (DL)",
   "virus"="virus (DL)",
   "fusion"="fusion (DL)"
)

revalue <- function(values, mappings){
   sapply(values, function(value){
      if(value %in% names(mappings))
         return(mappings[[value]])
      
      return(value)
   })
}

LOW_SNV_COUNT_THRES <- 50

PROBABILITIES_DECIMAL_PLACES <- 2

SIGNATURE_VALUES_DECIMAL_PLACES <- 0
SIGNATURE_PERC_SIGNIF_DIGITS <- 2

EVENT_FEATURES_DECIMAL_PLACES <- 1
EVENT_FEATURE_ABSENT_VALUE <- -0.00000001

CV_PERFORMANCE_DECIMAL_PLACES <- 2

NUMBER_FORMATTING_SIGNATURE_QUANTILES <- list(
   list(lower=10000,  upper=Inf,    func=function(x) formatC(x, format="e", digits=0)),
   list(lower=10,     upper=10000,  func=function(x) round(x, digits=0)),
   list(lower=1,      upper=10,     func=function(x) round(x, digits=1)),
   list(lower=-1,     upper=1,      func=function(x) round(x, digits=2)),
   list(lower=-10,    upper=-1,     func=function(x) round(x, digits=1)),
   list(lower=-10000, upper=-10,    func=function(x) round(x, digits=0)),
   list(lower=-Inf,   upper=-10000, func=function(x) formatC(x, format="e", digits=0))
)

NUMBER_FORMATTING_ODDS_RATIOS <- list(
   list(lower=10000, upper=Inf,   func=function(x) formatC(x, format="e", digits=1)),
   list(lower=10,    upper=10000, func=function(x) round(x, digits=0)),
   list(lower=1,     upper=10,    func=function(x) round(x, digits=1)),
   list(lower=0.01,  upper=1,     func=function(x) signif(x, digits=1)),
   list(lower=-Inf,  upper=0.01,  func=function(x) formatC(x, format="e", digits=0))
)

format_numbers <- function(values, configs){
   
   for(config in configs){
      if(!all(c("lower","upper","func") %in% names(config))){
         stop("`configs` must be list containing items in the format: `list(lower={numeric}, upper={numeric}, func={function})`")
      }
   }
   
   sapply(values, function(value){
      for(config in configs){
         if(value >= config$lower & value < config$upper){
            return(config$func(value))
         }
      }
      return(value)
   })
}

## Heatmap of probs ================================
get_plot_data_probs <- function(vis_data){

   plot_data <- subset(vis_data, data_type=="prob")

   plot_data$row_label <- revalue(plot_data$clf_name, MAPPINGS_CLASSIFIER_NAMES)
   plot_data$row_group <- toupper(plot_data$clf_group)

   ## Data label --------------------------------
   plot_data$data_label <- ifelse(
      !is.na(plot_data$data_value),
      round(plot_data$data_value, PROBABILITIES_DECIMAL_PLACES),
      ""
   )

   ## Remove rows for clfs without data --------------------------------
   clfs_missing_data <- sapply(
      split(plot_data$data_value, plot_data$clf_name),
      function(x) all(is.na(x))
   )
   clfs_missing_data <- names(clfs_missing_data)[clfs_missing_data]

   plot_data <- subset(plot_data, !(clf_name %in% clfs_missing_data))

   ## Get supertype probs --------------------------------
   plot_data$data_label_nudge_x <- 0

   split_plot_data_by_supertype <- function(plot_data){
      df <- subset(plot_data, grepl("combined", row_label, ignore.case=TRUE))
      split_values <- interaction(df$row_label, df$cancer_supertype)
      split(df, split_values)
   }

   plot_data_split <- split_plot_data_by_supertype(plot_data)

   plot_data_supertypes <- lapply(plot_data_split, function(df){
      df$row_label <- paste0("[GROUP] ", df$row_label)

      if(nrow(df)==1) return(df)

      prob_sum <- sum(df$data_value, na.rm=TRUE)
      prob_sum <- round(prob_sum, PROBABILITIES_DECIMAL_PLACES)

      ## Fill value
      df$data_value <- prob_sum

      ## Position label in the middle of row group
      middle_subclass <- ceiling(nrow(df)/2)
      df$data_label <- ""
      df$data_label[middle_subclass] <- prob_sum
      df$data_label_nudge_x <- if(nrow(df) %% 2 == 0){ 0.5 } else { 0 }

      return(df)
   })

   plot_data_supertypes <- do.call(rbind, unname(plot_data_supertypes))

   plot_data <- rbind(plot_data_supertypes, plot_data)

   ## Cell border color --------------------------------
   plot_data$cell_color <- ifelse(
      plot_data$row_label %in% plot_data_supertypes$row_label,
      NA,
      "darkgrey"
   )

   return(plot_data)
}

plot_heatmap <- function(
   plot_data,
   x_axis_position="top",
   y_strip_whitespaces=21,
   show_first_row_hline=FALSE,

   cell_color=NULL,
   data_label_nudge_x=NULL,

   simplify_cancer_type_labels=TRUE,
   cancer_type_label_wrap_width=20,

   y_labels_remap=waiver()
){
   ## Checks --------------------------------
   if(!(x_axis_position %in% c("top","bottom","none")))
      stop("`x_axis_position` must be 'top', 'bottom' or 'none'")

   ## Format labels --------------------------------
   ## y facet
   set_text_width <- function(x, whitespaces){

      white_space <- paste(rep(" ", whitespaces), collapse="")

      ## Use the above and below line to insert white spaces so that:
      ## - The width is constant
      ## - The target text is vertically centered
      x <- sprintf("%s\n%s\n%s", white_space, x, white_space)

      return(x)
   }

   if(!is.null(y_strip_whitespaces)){
      plot_data$row_group <- set_text_width(plot_data$row_group, y_strip_whitespaces)
   }

   ## Cancer type labels
   format_cancer_type_labels <- function(cancer_types){
      if(simplify_cancer_type_labels){
         cancer_types <- CANCER_TYPE_METADATA[cancer_types, "subtype"]
      }

      if(!is.null(cancer_type_label_wrap_width)){
         cancer_types <- str_wrap(
            as.character(cancer_types),
            whitespace_only=FALSE,
            width=cancer_type_label_wrap_width
         )
      }

      return(cancer_types)
   }

   ## Hide x facet strips where supertype only has itself as the subtype --------------------------------
   get_facet_strip_theme <- function(){

      subtypes_per_supertype <- CANCER_TYPE_METADATA[
         !duplicated(CANCER_TYPE_METADATA$supertype),
         "subtypes_per_supertype"
      ]

      supertype_has_subtypes <- subtypes_per_supertype > 1

      strip_fill <- "lightgrey"
      strip_color <- "black"

      strip_theme <- strip_themed(
         background_x = elem_list_rect(
            fill=ifelse(supertype_has_subtypes, strip_fill, "white"),
            color=ifelse(supertype_has_subtypes, strip_color, "#00000000")
         ),

         text_x = elem_list_text(
            color=ifelse(supertype_has_subtypes, "black", "#00000000")
         ),

         background_y=elem_list_rect(
            fill=strip_fill,
            color=strip_color
         )
      )

      return(strip_theme)
   }

   facet_strip_theme <- get_facet_strip_theme()

   ## Plot --------------------------------
   ## Convert character vectors to factors
   for(colname in colnames(plot_data)){
      if(is.character(plot_data[[colname]])){
         plot_data[[colname]] <- factor(plot_data[[colname]], unique(plot_data[[colname]]))
      }
   }

   p <- ggplot(plot_data, aes(y=row_label, x=cancer_type)) +
      facet_grid2(
         row_group ~ cancer_supertype,
         scales="free", space="free",
         strip=if(x_axis_position=="none") strip_vanilla() else facet_strip_theme,
         switch=if(x_axis_position=="bottom") "both" else "y"
      ) +

      geom_tile(
          aes(fill=data_value),
          color = if(is.null(cell_color)) NA else as.character(cell_color)
      ) +

      geom_text(
          aes(label=data_label), size=2.7,
          nudge_x = if(is.null(data_label_nudge_x)) 0 else plot_data[[data_label_nudge_x]]
      ) +

      coord_cartesian(expand=FALSE) +
      scale_y_discrete(limits=rev, labels=y_labels_remap) +
      scale_x_discrete(
         position=if(x_axis_position=="bottom") "bottom" else "top",
         name="Cancer supertype (strips) and subtype (labels)",
         labels=format_cancer_type_labels
      ) +

      theme_bw() +
      theme(
         panel.grid=element_blank(),
         panel.spacing=unit(3, "pt"),
         panel.border=element_rect(linewidth=0.75),
         #strip.placement.x="outside", ## patchwork::wrap_plots() fails with this
         strip.placement="outside",
         strip.text.y.right=element_text(angle=0),
         strip.text.y.left=element_text(angle=0, hjust=1),
         axis.text.x.top=element_text(angle=90, hjust=0, vjust=0.5),
         axis.text.x.bottom=element_text(angle=90, hjust=1, vjust=0.5),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         legend.margin=margin(c(0,0,0,0)),
         legend.justification=c("left","top"),
         plot.title=element_text(face="bold"),
         plot.margin = margin(0,0,20,0, "pt")
      )

   ## Optional plot modifications --------------------------------
   ## Add horizontal line on first row
   get_plot_data_first_row_hline <- function(plot_data){

      plot_data_hline <- plot_data[c("row_label", "row_group", "cancer_type", "cancer_supertype")]
      plot_data_hline <- plot_data_hline[!duplicated(interaction(plot_data_hline$cancer_type, plot_data_hline$row_label)),]

      facet_n_rows <- unique(plot_data_hline[c("row_label", "row_group")])
      facet_n_rows <- unclass(table(facet_n_rows$row_group))

      plot_data_hline$yintercept <- facet_n_rows[ plot_data_hline$row_group ]

      return(plot_data_hline)
   }


   if(show_first_row_hline){
      plot_data_hline <- get_plot_data_first_row_hline(plot_data)
      p <- p + geom_hline(data=plot_data_hline, mapping=aes(yintercept=yintercept-0.5), linetype="A1")
   }

   ## Remove x axis
   if(x_axis_position=="none"){
      p <- p +
         theme(
            axis.text.x.bottom=element_blank(),
            axis.ticks.x.bottom=element_blank(),
            axis.text.x.top=element_blank(),
            axis.ticks.x.top=element_blank(),
            strip.background.x = element_blank(),
            strip.text.x = element_blank()
         )
   }

   ## Store heat map dimensions within plot --------------------------------
   p$ncols <- length(unique(plot_data$cancer_type))
   p$nrows <- length(unique(plot_data$row_label))

   return(p)
}

get_snv_count <- function(VIS_DATA){
   snv_count <- subset(VIS_DATA,  str_ends(feat_name, "snv_count"), feat_value)[1,1]
   if(is.na(snv_count)){
      stop("`snv_count` was not found as a feature in `VIS_DATA`")
   }
   return(snv_count)
}

plot_probs <- function(vis_data){
   plot_data <- get_plot_data_probs(vis_data)

   snv_count <- get_snv_count(vis_data)
   if(snv_count < LOW_SNV_COUNT_THRES){
      warning_string <- sprintf(" | WARNING: Sample has <%s SNVs. Predictions may be unreliable", LOW_SNV_COUNT_THRES)
   } else {
      warning_string <- ""
   }

   plot_heatmap(
      plot_data,
      show_first_row_hline = TRUE,
      cell_color = plot_data$cell_color,
      data_label_nudge_x = "data_label_nudge_x"
   ) +
   scale_fill_gradient(
      low="white", high="mediumseagreen", na.value="grey95", limits=c(0,1),
      guide=guide_colorbar(frame.colour="black", ticks.colour="black", barheight=3.5, frame.linewidth=0.25)
   ) +
   labs(
      title=paste0("Probabilities by classifier", warning_string),
      subtitle="Cancer group (strips) and subtype (label)",
      fill="Probability",
   )
}


## Plotting other data_types ================================
plot_cv_performance <- function(vis_data){
   plot_data <- subset(vis_data, data_type=="cv_performance")

   plot_data$row_group <- "training set"
   plot_data$row_label <- plot_data$feat_name

   ## Data label --------------------------------
   is_n_samples_row <- plot_data$feat_name=="n_total"

   plot_data$data_label <- with(plot_data, {
      data_label <- as.character(round(data_value, CV_PERFORMANCE_DECIMAL_PLACES))
      data_label[is_n_samples_row] <- data_value[is_n_samples_row]
      return(data_label)
   })

   ## Scale values for n_samples to max of 1 --------------------------------
   plot_data$data_value <- with(plot_data, {
      n_samples <- data_value[is_n_samples_row]
      n_samples_scaled <- n_samples / max(n_samples)

      data_value[is_n_samples_row] <- n_samples_scaled

      return(data_value)
   })

   ## Plot --------------------------------
   y_labels_remap <- c(
      n_total="Total no. of samples",
      recall="Recall (prop. of total correct)",
      precision="Precision (prop. of predicted correct)"
   )

   plot_heatmap(
      plot_data,
      show_first_row_hline = TRUE,
      x_axis_position="none",
      cell_color="lightgrey",
      y_labels_remap = y_labels_remap

   ) +
   scale_fill_gradientn(
      colors=c("white", "#fffedc", "#76c99b"), na.value="grey95", limits=c(0,1),
      guide=guide_colorbar(frame.colour="black", ticks.colour="black", barheight=3.5, frame.linewidth=0.25)
   ) +
   labs(
      title=sprintf("%s: Training set performance", MAPPINGS_CLASSIFIER_NAMES[["dna_combined"]]), 
      fill="Metric value"
   ) +
   theme(legend.justification=c("left", "bottom"))
}

plot_signatures <- function(vis_data){

   plot_data <- subset(vis_data, data_type=="sig_quantile")

   plot_data$row_group <- "signature"
   
   plot_data$data_label <- format_numbers(plot_data$data_value, NUMBER_FORMATTING_SIGNATURE_QUANTILES)

   ## Make row labels --------------------------------
   snv_count <- get_snv_count(vis_data)
   
   plot_data$row_label <- with(plot_data, {
      perc <- round((feat_value / snv_count) * 100, SIGNATURE_PERC_SIGNIF_DIGITS)
      feat_value <- round(feat_value, SIGNATURE_VALUES_DECIMAL_PLACES)
      paste0(gsub("_", " ", feat_name), " = ", feat_value, " (", perc, "%)")
   })

   ## Discretize quantiles --------------------------------
   quantile_info <- list(
      breaks=c(-Inf, -0.004, 0.004, 0.95, 1.2, Inf),
      labels=c("<0 (below expected range)", "0", ">0 - 0.95 (in expected range)", "0.95 - 1.2 (above expected range)", ">1.2 (well above expected range)"),
      colors=c("lightcyan2", "grey95", "#C9E7CD", "#f5dfdf","#e8b6b6") ## white, green, light green, red
   )

   plot_data$data_value <- cut(
      plot_data$data_value,
      breaks=quantile_info$breaks,
      labels=quantile_info$labels
   )

   ## Plot --------------------------------
   plot_heatmap(plot_data, x_axis_position="none", cell_color="grey50") +
   scale_fill_manual(
      values=structure(
         quantile_info$colors[quantile_info$colors != "grey95"], 
         names=quantile_info$labels[quantile_info$labels != "0"]
      ),
      guide=guide_legend(override.aes=list(colour="black", linetype=1, linewidth=0.25)),
      na.value="grey95", drop=TRUE
   ) +
   labs(
      fill="Quantile in subtype cohort", 
      title=sprintf("%s: Mutational signatures", MAPPINGS_CLASSIFIER_NAMES[["snv96"]])
   )
}

plot_feat_contrib <- function(vis_data){

   plot_data <- subset(vis_data, data_type=="feat_contrib")
   
   ## Format labels --------------------------------
   ## Parse feature names
   affixes <- str_split(plot_data$feat_name, "[.]", n=2, simplify=TRUE)
   feat_types <- affixes[,1]
   feat_names <- affixes[,2]
   rm(affixes)

   ## Format odds ratios
   odds_ratios <- format_numbers(plot_data$data_value, NUMBER_FORMATTING_ODDS_RATIOS)
   odds_ratios[odds_ratios==1] <- "" ## odds of 1 == log(odds) of 0 -> unimportant, thus hide
   odds_ratios <- sub("[-]0", "-", odds_ratios) ## Remove leading 0 from negative exponent
   odds_ratios <- sub("[+]0", "", odds_ratios) ## Remove leading 0 and sign from positive exponent
   plot_data$data_label <- odds_ratios

   ## Make feature names more human readable
   feat_names <- sapply(1:nrow(plot_data), function(i){
      feat_type <- feat_types[i]
      feat_name <- feat_names[i]
      
      if(feat_type %in% c("driver"))
         return(sub(".", " ", feat_name, fixed=TRUE))
      
      if(feat_type == "fusion")
         return(sub("_", "-", feat_name, fixed=TRUE))
      
      if(feat_type %in% c("trait","tmb","sv"))
         return(gsub("_", " ", feat_name, fixed=TRUE))
      
      return(feat_name)
   })
   
   ## Format feature values
   raw_feat_values <- plot_data$feat_value
   feat_values <- round(raw_feat_values, digits=EVENT_FEATURES_DECIMAL_PLACES)
   feat_values[raw_feat_values==0] <- "0" ## Use integer representation
   feat_values[raw_feat_values==1] <- "1" ## Use integer representation
   feat_values[raw_feat_values==EVENT_FEATURE_ABSENT_VALUE] <- "absent"
   feat_values <- ifelse(feat_types=="trait", feat_values==1, feat_values) ## Convert to boolean
   
   plot_data$row_label <- ifelse(
      feat_types=="driver",
      sprintf("%s (%s)", feat_names, feat_values),
      sprintf("%s = %s", feat_names, feat_values)
   )

   plot_data$row_group <- revalue(feat_types, MAPPING_EVENT_NAMES)
   
   ## Legend breaks --------------------------------
   exponent_range <- log10(max(plot_data$data_value)) - log10(min(plot_data$data_value))
   legend_breaks <- 10^seq(from=-10, to=10, by=if(exponent_range < 5) 1 else 2)

   ## Plot --------------------------------
   plot_heatmap(
      plot_data,
      x_axis_position="none",
      cell_color="lightgrey"
   ) +
   scale_fill_gradient2(
      low="indianred", high="mediumseagreen", trans="log",
      breaks=legend_breaks,
      labels = function(x) sprintf("%g", x),
      guide=guide_colorbar(frame.colour="black", ticks.colour="black", barheight=5, frame.linewidth=0.25)
   ) +
   labs(
      fill="Odds ratio\n(subtype / not subtype)", 
      y="Feature", 
      title=sprintf("%s: Feature contributions", MAPPINGS_CLASSIFIER_NAMES[["event"]])
   )
}


## Combine and export ================================
plot_one_sample <- function(VIS_DATA, sample_id){
   
   vis_data <- VIS_DATA[VIS_DATA$sample_id==sample_id | VIS_DATA$data_type=="cv_performance",]
   
   ## Plots --------------------------------
   plots <- list()
   
   plots$probs <- plot_probs(vis_data)
   plots$signatures <- plot_signatures(vis_data)
   plots$feat_contrib <- plot_feat_contrib(vis_data)
   
   if(!("cv_performance" %in% vis_data$data_type)){
      stop("CV performance data is missing from vis data")
   }
   plots$cv_performance <- plot_cv_performance(vis_data)
   
   ## Combine --------------------------------
   ## Heights
   heights <- sapply(plots, function(p){ p$nrows })
   heights <- heights[heights>0] ## Remove 0 height plots
   
   ## Widths
   widths <-  sapply(plots, function(p){ p$ncols })
   
   plots_combined <- patchwork::wrap_plots(plots, ncol=1, heights=heights)
   
   plots_combined$heights <- heights
   plots_combined$widths <- widths
   
   return(plots_combined)
}

write_plots <- function(VIS_DATA){
   
   sample_ids <- unique(VIS_DATA[VIS_DATA$data_type=="prob","sample_id"])
   
   ## Calculate dimensions --------------------------------
   if(length(sample_ids)==1){
      
      p <- plot_one_sample(VIS_DATA, sample_ids[1])
      
      ## Fixed height/width account for axis labels
      fixed_height <- 3
      fixed_width <- 6
      
      height_scaling <- 0.25
      width_scaling <- 0.35
      
      height <- sum(p$heights) * height_scaling + fixed_height
      width  <- max(p$widths) * width_scaling + fixed_width
      
   } else {
      ## Use constant dimensions for multi sample plotting
      height <- 12.5
      width <- 20
   }
   
   ## Open output file connection --------------------------------
   if(grepl(".pdf$", PLOT_PATH)){
      pdf(PLOT_PATH, width, height)
   } else if(grepl(".png$", PLOT_PATH)){
      if(length(sample_ids)>1) stop("Plot output file extension must be .pdf for multi-sample plotting")
      png(PLOT_PATH, width, height, units="in", res=300)
   } else {
      stop("Plot output file must have .pdf or .png as the extension")
   }
   
   ## Write to file --------------------------------
   if(length(sample_ids)==1){
      plot(p)
   } else {
      for(sample_id in sample_ids){
         message("Plotting sample: ", sample_id)
         p <- 
            plot_one_sample(VIS_DATA, sample_id) + 
            plot_annotation(
               paste0("Sample: ", sample_id),
               theme = theme(plot.title=element_text(face="bold"))
            )
         plot(p)
      }
   }
   
   invisible(dev.off())
}

write_plots(VIS_DATA)
