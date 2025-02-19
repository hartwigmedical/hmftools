library(dplyr)
library(tidyr)
library(ggplot2)

# Parse and check inputs
args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) < 3)
{
  print("Requires arguments 1=Directory, 2=SampleId, 3=MsStatsTableFile")
  stop()
}

outDir <- args[1]
sampleId <- args[2]
msStatsTableFile <- args[3]

if (!file.exists(msStatsTableFile))
{
  print(sprintf('Missing microsatellite stats sample data file: %s', msStatsTableFile))
  stop()
}

df <- read.table(msStatsTableFile, header = TRUE, sep = "\t", check.names = FALSE) %>%
  rename("reads" = "readCount")

# fill in missing rows
df <- df %>%
  merge(expand.grid(unit = unique(df$unit), consensusType = unique(df$consensusType), numUnits = c(4:20)), all=TRUE) %>%
  replace(is.na(.), 0)

# https://stackoverflow.com/questions/58515587/heat-map-with-additional-values-in-r
# change to long form
df <- df %>% 
  pivot_longer(
    cols = starts_with('count'), 
    names_to = "jitter",
    values_to = "count")

df$jitter <- as.character(gsub("count", "", df$jitter))
df$rate <- pmax(df$count / df$reads * 100, 0.001)

df$rate <- replace(df$rate, is.na(df$rate), 0)

# make sure we have a row for read count for every unit, consensusType, numUnits tuple
df <- df %>%
  group_by(unit, consensusType, numUnits) %>%
  summarise(jitter="reads", rate=first(reads), .groups = 'keep') %>%
  rbind(df)

# need to give it the sorting of the jitter values
df$jitter <- factor(df$jitter, levels=c(sprintf("%+d", c(-10:10)), "reads"))

units <- c("A/T", "C/G", "AT/TA", "AG/GA/CT/TC", "AC/CA/GT/TG", "CG/GC", "3bp repeat", "4bp repeat", "5bp repeat")

# Create a list to store individual plots
plot_list <- list()

format_read_counts <- function(read_counts, scientific_cutoff=10000){
  ifelse(
    read_counts < scientific_cutoff,
    format(read_counts, big.mark=",", scientific=FALSE),
    sub("E[+]0", "E", sprintf("%.2G", read_counts))
  )
}

consensusTypes <- unique(df$consensusType)
for(msConsensusType in consensusTypes)
{
  msStatsHeatmap <- file.path(outDir, sprintf('%s.%s.ms_stats.png', sampleId, tolower(msConsensusType)))
  png(file = msStatsHeatmap, res = 140, height = 2200, width = 4000)

  for (msUnit in units) {
    plot_key <- sprintf('unit(%s) consensusType(%s)', msUnit, msConsensusType)
    print(plot_key)
    jitter_df <- df %>%
      filter(unit == msUnit & consensusType == msConsensusType, numUnits %in% 4:20) %>%
      select(-c(unit))
    
    # Create a plot for each unit
    p <- ggplot(jitter_df, aes(x = jitter, y = reorder(numUnits, desc(numUnits)))) +
      geom_tile(aes(fill = rate), color = "white") +
      scale_fill_gradient(low = "white", high = "deepskyblue", trans = "log", limits=c(0.001,100), guide="none", na.value="lightgrey") +
      geom_text(aes(label = sprintf("%.2f", rate)), size=3, subset(jitter_df, jitter != "reads" & reads > 0)) +
      
      geom_text(aes(label = format_read_counts(rate)), size=3, subset(jitter_df, jitter == "reads")) +
      
      theme_minimal() +
      labs(title = sprintf('unit = %s consensusType = %s', msUnit, msConsensusType), x="jitter", y="num units")
    
    # Add the plot to the list
    plot_list[[msUnit]] <- p
  }

  # Arrange the plots in a 3x3 grid
  gridExtra::grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)
  dev.off()
}