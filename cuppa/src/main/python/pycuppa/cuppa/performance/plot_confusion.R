options(max.print=500)
options(stringsAsFactors=FALSE)

library(ggplot2)
library(reshape2)
suppressWarnings(library(patchwork))

## Args ================================
args <- commandArgs(trailingOnly = TRUE)

PRED_SUMM_PATH <- args[1]
CLF_NAME <- args[2]
PLOT_PATH <- args[3]
PLOT_WIDTH <- args[4]
PLOT_HEIGHT <- args[5]

if(FALSE){
   PRED_SUMM_PATH <- "pred_summ.tsv.gz"
   CLF_NAME <- "dna_combined"
   PLOT_PATH <- "plot.pdf"
   PLOT_WIDTH <- 14
   PLOT_HEIGHT <- 10
}

if(is.na(PLOT_PATH)){
   stop("PLOT_PATH must be provided")
}

## Load data ================================
pred_summ <- read.delim(PRED_SUMM_PATH)

if(!("actual_class" %in% colnames(pred_summ))){
   stop("Cannot plot confusion matrix `actual_class` column is absent")
}

if(!(CLF_NAME %in% pred_summ$clf_name)){
   stop(sprintf("'%s' is not a valid `clf_name`", CLF_NAME))
}

## Transform data ================================
get_confusion <- function(pred_summ){
   df <- subset(pred_summ, clf_name==CLF_NAME)
   
   unique_classes <- sort(unique(c(".All", df$pred_class_1, df$actual_class)))
   
   df$pred_class_1 <- factor(df$pred_class_1, unique_classes)
   df$actual_class <- factor(df$actual_class, unique_classes)
   
   confusion <- unclass(table(
      pred_class = df$pred_class_1, 
      actual_class = df$actual_class
   ))
   
   confusion[,".All"] <- NA
   confusion[".All",] <- NA
   
   return(confusion)
}

get_performance <- function(confusion){
   
   stats <- data.frame(
      n_total = colSums(confusion, na.rm = TRUE),
      n_correct = diag(confusion),
      n_predicted = rowSums(confusion, na.rm = TRUE)
   )
   
   stats[".All",] <- colSums(stats)
   stats$n_correct[1] <- sum(stats$n_correct, na.rm = TRUE)
   stats$n_predicted[1] <- NA
   
   stats$recall <- stats$n_correct / stats$n_total
   stats$precision <- stats$n_correct / stats$n_predicted
   
   return(stats)
}

confusion <- get_confusion(pred_summ)
performance <- get_performance(confusion)

## Plotting ================================
FILL_COLORS <- c('#225EA8', '#1D91C0', '#41B6C4', '#7FCDBB', '#C7E9B4', '#EDF8B1', '#FFFFD9')

plot_heatmap <- function(
   plot_data, 
   x, y, fill, label,
   direction="horizontal",
   x_title = "Actual class", y_title = "Predicted class",
   hide_x = FALSE, hide_y = FALSE
){
   if(FALSE){
      x="actual_class"
      y="pred_class"
      fill="prop"
      label="count"
      
      x="class"
      y="metric"
      fill="fill_value"
      label="label_value"
      
      direction = "vertical"
      x_title = "Actual class"
      y_title = "Predicted class"
   }
   
   if(direction == "horizontal"){
      colname_mappings <- c(x=x, y=y, fill=fill, label=label)
   } else if (direction == "vertical"){
      colname_mappings <- c(x=y, y=x, fill=fill, label=label)
   } else {
      stop("`direction` must be 'horizontal' or 'vertical'")
   }
   
   for(dest_name in names(colname_mappings)){
      orig_name <- colname_mappings[[dest_name]]
      colnames(plot_data)[colnames(plot_data)==orig_name] <- dest_name
   }
   
   p <- ggplot(plot_data, aes(y=y, x=x)) +
      
      geom_tile(aes(fill=fill), color="black") +
      geom_text(aes(label=label), size=2.7) +
      scale_fill_gradientn(
         name = "Proportion",
         colors=FILL_COLORS, limits=c(0,1),
         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
      ) +
      
      coord_cartesian(expand = FALSE) +
      scale_x_discrete(name=x_title, position = "top") +
      scale_y_discrete(name=y_title, limits = rev) +
      
      theme_bw() +
      theme(
         panel.grid = element_blank(),
         axis.text.x.bottom = element_text(angle=90, hjust=1, vjust=0.5),
         axis.text.x.top = element_text(angle=90, hjust=0, vjust=0.5),
         plot.margin = margin(0,0,0,0,"pt")
      )
   
   if(hide_x){
      p <- p + theme(
         axis.title.x = element_blank(),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.title.x.top = element_blank(),
         axis.text.x.top = element_blank(),
         axis.ticks.x.top = element_blank()
      )
   }
   
   if(hide_y){
      p <- p + theme(
         axis.title.y = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks.y = element_blank()
      )
   }
   
   p$height <- length(unique(plot_data$y))
   p$width <- length(unique(plot_data$x))
   
   return(p)
}

plot_confusion <- function(confusion){
   count <- melt(confusion, value.name = "count")
   prop <- melt(confusion / colSums(confusion, na.rm = TRUE), value.name = "prop")
   
   plot_data <- count
   plot_data$prop <- prop$prop
   
   performance <- get_performance(confusion)
   pancancer_row <- plot_data$pred_class == ".All" & plot_data$actual_class == ".All"
   plot_data[pancancer_row, "count"] <- performance[".All", "n_correct"]
   plot_data[pancancer_row, "prop"] <- performance[".All", "recall"]
  
   plot_heatmap(plot_data,  x="actual_class", y="pred_class", fill="prop", label="count")
}

plot_stats <- function(performance, metrics=NULL, direction="horizontal", hide_x=FALSE, hide_y=FALSE){
   
   if(FALSE){
      metrics <- c("recall","n_total")
      direction="horizontal"
      hide_x=FALSE
      hide_y=FALSE
   }
   
   plot_data <- melt(as.matrix(performance))
   colnames(plot_data) <- c("class", "metric", "value")
   
   maxes <- apply(performance, 2, max, na.rm=TRUE)
   plot_data$max <- maxes[as.character(plot_data$metric)]
   
   plot_data$is_count <- plot_data$metric %in% c("n_total","n_correct","n_predicted")
   
   plot_data$fill_value <- ifelse(
      plot_data$is_count,
      log10(plot_data$value) / log10(plot_data$max),
      plot_data$value / plot_data$max
   )
   
   ## Labels
   plot_data$label_value <- ifelse(
      plot_data$metric %in% c("recall","precision"),
      round(plot_data$value, 2), 
      plot_data$value
   )
   
   ##
   if(!is.null(metrics)){
      plot_data <- plot_data[plot_data$metric %in% metrics, ]
      plot_data$metric <- factor(plot_data$metric, metrics)
   }
   
   plot_heatmap(
      plot_data, 
      x="class", y="metric", fill="fill_value", label="label_value", 
      direction=direction,
      x_title = "",
      y_title = "",
      hide_x=hide_x, 
      hide_y=hide_y
   )
}

## Export ================================
combine_plots <- function(confusion, performance){
   plots <- list()
   
   plots$confusion <- plot_confusion(confusion)
   plots$row_stats <- plot_stats(performance, metrics=c("precision", "n_correct", "n_predicted"), direction = "vertical", hide_y=TRUE)
   plots$col_stats <- plot_stats(performance, metrics=c("recall", "n_correct", "n_total"), direction = "horizontal", hide_x=TRUE)
   plots$spacer <- plot_spacer()
   
   widths <- c(plots$confusion$width, plots$row_stats$width)
   heights <- c(plots$confusion$height, plots$col_stats$height)
   
   p <- wrap_plots(plots, guides = "collect", widths=widths, heights=heights)
   p$width <- sum(widths)
   p$height <- sum(heights)
   
   return(p)
}

p <- combine_plots(confusion, performance)
width <- if(is.na(PLOT_WIDTH)) p$width*0.35 else PLOT_WIDTH
height <- if(is.na(PLOT_HEIGHT)) p$height*0.25 else PLOT_HEIGHT

suppressWarnings({
   pdf(PLOT_PATH, width, height)
   plot(p)
   dev.off()
})



