options(max.print=500)

library(ggplot2)
library(ggh4x)
library(stringr)
library(patchwork)

## Args ================================
args <- commandArgs(trailingOnly = TRUE)

VIS_DATA_PATH <- args[1]
PLOT_PATH <- args[2]

if(FALSE){
   VIS_DATA_PATH <- "~/Desktop/cuppa_vis_data.tsv"
   PLOT_PATH <- "~/Desktop/cuppa_vis.png"
}

## Load data ================================
VIS_DATA <- read.delim(VIS_DATA_PATH)

as_factor_unsorted <- function(x){
   if(is.vector(x)){
      return(factor(x, unique(x)))
   }

   if(is.data.frame(x)){
      for(colname in colnames(x)){
         if(is.character(x[[colname]])){
            x[[colname]] <- factor(x[[colname]], unique(x[[colname]]))
         }
      }
      return(x)
   }

   stop("`x` must be a vector or dataframe")
}

get_cancer_type_metadata <- function(VIS_DATA){

   cancer_types <- unique(VIS_DATA$cancer_type)

   prefix <- str_extract(cancer_types, "^[A-Za-z/ ]+")
   suffix <- str_replace(cancer_types, "^[A-Za-z/ ]+: ", "")

   df <- data.frame(
      supertype = prefix,
      subtype = suffix,

      row.names = cancer_types
   )
   df <- as_factor_unsorted(df)

   subtypes_per_supertype <- sapply(split(df$subtype, df$supertype), length)
   df$subtypes_per_supertype <- subtypes_per_supertype[df$supertype]

   return(df)
}

CANCER_TYPE_METADATA <- get_cancer_type_metadata(VIS_DATA)

VIS_DATA$cancer_supertype <- CANCER_TYPE_METADATA[VIS_DATA$cancer_type, "supertype"]


## Heatmap of probs ================================
get_plot_data_probs <- function(VIS_DATA, data_value_rounding = 2){

   plot_data <- subset(VIS_DATA, data_type=="prob")

   plot_data$row_label <- toupper(plot_data$clf_name)
   plot_data$row_group <- toupper(plot_data$clf_group)

   ## Data label --------------------------------
   plot_data$data_label <- ifelse(
      !is.na(plot_data$data_value),
      round(plot_data$data_value, data_value_rounding),
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
      #df = plot_data_split[["COMBINED.Bone/Soft tissue"]]
      df$row_label <- paste0("[GROUP] ", df$row_label)

      if(nrow(df)==1) return(df)

      prob_sum <- sum(df$data_value, na.rm=TRUE)
      prob_sum <- round(prob_sum, 2)

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

   plot_data <- rbind(
      plot_data_supertypes,
      plot_data
   )

   ## Cell border color --------------------------------
   ## Remove color for
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

   if(FALSE){
      plot_data <- get_plot_data_probs(VIS_DATA)
      x_axis_position = "top"
      y_strip_whitespaces = NULL
      show_first_row_hline = TRUE

      cell_color = plot_data$cell_color
      data_label_nudge_x = "data_label_nudge_x"

      simplify_cancer_type_labels=TRUE
      cancer_type_label_wrap_width=20
   }

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
      x <- sprintf(
         "%s\n%s\n%s",
         white_space, x, white_space
      )

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
   plot_data <- as_factor_unsorted(plot_data)

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

      # scale_fill_gradient(
      #    low="white", high="mediumseagreen", na.value="grey95", limits=c(0,1),
      #    guide=guide_colorbar(frame.colour="black", ticks.colour="black", barheight=5)
      # ) +

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

plot_probs <- function(VIS_DATA){
   plot_data <- get_plot_data_probs(VIS_DATA)

   plot_heatmap(
      plot_data,
      show_first_row_hline = TRUE,
      cell_color = plot_data$cell_color,
      data_label_nudge_x = "data_label_nudge_x"
   ) +
   scale_fill_gradient(
      low="white", high="mediumseagreen", na.value="grey95", limits=c(0,1),
      guide=guide_colorbar(frame.colour="black", ticks.colour="black", barheight=3.5)
   ) +
   labs(
      title="Probabilities by classifier",
      subtitle="Cancer group (strips) and subtype (label)",
      fill="Probability",
   )
}


## Plotting other data_types ================================
plot_cv_performance <- function(VIS_DATA, data_value_rounding = 2){
   plot_data <- subset(VIS_DATA, data_type=="cv_performance")

   plot_data$row_group <- "training_set"
   plot_data$row_label <- plot_data$feat_name

   ## Data label --------------------------------
   is_n_samples_row <- plot_data$feat_name=="n_samples"

   plot_data$data_label <- with(plot_data, {
      data_label <- as.character(round(data_value, 2))
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
         n_samples="Total no. of samples",
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
      guide=guide_colorbar(frame.colour="black", ticks.colour="black", barheight=3.5)
   ) +
   #guides(fill=FALSE) +
   labs(title="DNA_COMBINED: Training set performance", fill="Metric value") +
   theme(legend.justification=c("left", "bottom"))
}


plot_signatures <- function(
   VIS_DATA,
   data_value_rounding = 2,
   feat_value_rounding = 0,
   feat_perc_signif = 1
){
   if(FALSE){
      data_value_rounding = 2
      feat_value_rounding = 0
      feat_perc_signif = 2
   }

   plot_data <- subset(VIS_DATA, data_type=="sig_quantile")

   plot_data$row_group <- "signature"
   plot_data$data_label <- round(plot_data$data_value, data_value_rounding)

   ## Make row labels --------------------------------
   snv_count <- subset(VIS_DATA,  str_ends(feat_name, "snv_count"), feat_value)[1,1]

   plot_data$row_label <- with(plot_data, {

      perc <- round(
         (feat_value / snv_count) * 100,
         feat_perc_signif
      )

      feat_value <- round(feat_value, feat_value_rounding)

      paste0(feat_name, " = ", feat_value, " (", perc, "%)")
   })

   ## Discretize quantiles --------------------------------
   quantile_info <- list(
      breaks=c(0, 0.95, 1.2, Inf),
      labels=c(">0 - 0.95 (in expected range)", "0.95 - 1.2 (above expected range)", ">1.2 (well above expected range)"),
      colors=c("#C9E7CD", "#f5dfdf","#e8b6b6") ## white, green, light green, red
   )

   plot_data$data_value <- cut(
      plot_data$data_value,
      breaks=quantile_info$breaks,
      labels=quantile_info$labels
   )

   ## Plot --------------------------------
   plot_heatmap(
      plot_data,
      x_axis_position="none",
      cell_color="grey50"
   ) +
   scale_fill_manual(
      values=structure(quantile_info$colors, names=quantile_info$labels),
      guide=guide_legend(override.aes=list(colour="black", linetype=1)),
      na.value="grey95", drop=FALSE
   ) +
   labs(fill="Quantile in subtype cohort", title="SNV96: Mutational signatures")
}

plot_feat_contrib <- function(
   VIS_DATA,
   feat_value_rounding = 2,
   small_data_value_rounding = 1,
   large_data_value_rounding = 0,
   large_data_value_thres = 100
){
   if(FALSE){
      feat_value_rounding = 2
      small_data_value_rounding = 1
      large_data_value_rounding = 0
      large_data_value_thres = 100
   }

   plot_data <- subset(VIS_DATA, data_type=="feat_contrib")

   ## Parse feature names --------------------------------
   affixes <- as.data.frame(do.call(
      rbind,
      str_split(plot_data$feat_name, '[.]', n=2)
   ))
   colnames(affixes) <- c("feat_type", "feat_basename")

   ## Format labels --------------------------------
   plot_data$row_group <- affixes$feat_type

   ## Odds
   plot_data$data_label <-  with(plot_data, {
      data_label <- round(data_value, small_data_value_rounding)
      data_label[data_label==1] <- ""

      is_large_value <- data_value >= large_data_value_thres
      data_label[is_large_value] <- round(data_value, large_data_value_rounding)[is_large_value]

      return(data_label)
   })

   ## Row values
   plot_data$row_label <- with(plot_data, {

      feat_value_label <- round(feat_value, feat_value_rounding)
      feat_value_label[feat_value==0] <- "0"
      feat_value_label[feat_value==1] <- "1"

      paste0(affixes$feat_basename, " = ", feat_value_label)
   })

   ## Plot --------------------------------
   plot_heatmap(
      plot_data,
      x_axis_position="none",
      cell_color="lightgrey"
   ) +
   scale_fill_gradient2(
      low="indianred", high="mediumseagreen", trans="log",
      breaks=10^(-4:4),
      labels = function(x) sprintf("%g", x),
      guide=guide_colorbar(frame.colour="black", ticks.colour="black", barheight=5)
   ) +
   labs(fill="Odds ratio\n(subtype / not subtype)", y="Feature", title="EVENT: Feature contributions")
}


## Combine plots ================================
cuppa_vis <- function(VIS_DATA, plot_path=NULL){

   ## Plots --------------------------------
   plots <- list()

   plots$probs <- plot_probs(VIS_DATA)
   plots$signatures <- plot_signatures(VIS_DATA)
   plots$feat_contrib <- plot_feat_contrib(VIS_DATA)

   if("cv_performance" %in% VIS_DATA$data_type){
      plots$cv_performance <- plot_cv_performance(VIS_DATA)
   } else {
      warning("CV performance data is missing from vis data and will not be plotted")
   }

   ## Combine --------------------------------
   ## Heights
   heights <- sapply(plots, function(p){ p$nrows })
   heights <- heights[heights>0] ## Remove 0 height plots

   ## Widths
   widths <-  sapply(plots, function(p){ p$ncols })

   plots_combined <- patchwork::wrap_plots(plots, ncol=1, heights=heights)

   ## Export --------------------------------
   if(is.null(plot_path)) return(plots_combined)

   fixed_height <- 3
   ggsave(
      plot=plots_combined,
      filename=plot_path,
      height=sum(heights) * 0.25 + fixed_height,
      width=max(widths) * 0.5,
      units="in"
   )
}

cuppa_vis(VIS_DATA, PLOT_PATH)
