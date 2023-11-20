library(dplyr)
library(tidyr)
library(ggplot2)
library(gtable)

# Parse and check inputs
args <- commandArgs(trailingOnly = TRUE)
print(args)

if (length(args) < 5)
{
  print("Requires arguments 1=PatientId 2=SampleId 3=Summary file 4=CN Segment file 5=Sample data output directory")
  stop()
}

patientId = args[1]
sampleId = args[2]
summaryFile = args[3]
cnSegmentFile = args[4]
outputDir <- args[5]

print(sprintf('Producing CN-GC Ratio plot for patient(%s) sample(%s), writing to output(%s)', patientId, sampleId, outputDir))

if (!file.exists(summaryFile))
{
  print(sprintf('Missing sample summary file: %s', summaryFile))
  stop()
}

if (!file.exists(cnSegmentFile))
{
  print(sprintf('Missing sample CN-GcRatio file: %s', cnSegmentFile))
  stop()
}

sampleSummary = read.csv(summaryFile,sep='\t')
cnGcRatioSegments = read.csv(cnSegmentFile,sep='\t')

if('SampleId' %in% colnames(sampleSummary))
{
    sampleSummary = sampleSummary %>% filter(SampleId==sampleId)
    cnGcRatioSegments = cnGcRatioSegments %>% filter(SampleId==sampleId)
}

if(nrow(sampleSummary) == 0 | nrow(cnGcRatioSegments) == 0)
{
  print(sprintf('Empty filtered input files'))
  stop()
}

cnFitIntercept=sampleSummary$CnFitIntercept
cnFitCoeff=sampleSummary$CnFitCoeff

# cnGcRatioSegments = merge(cnGcRatioSegments,sampleSummary %>% select(SampleId,CnFitIntercept,CnFitCoeff),by='SampleId',all.x=T)

cnGcRatioSegments = cnGcRatioSegments %>%
    mutate(GcRatioMedianFit=cnFitIntercept+cnFitCoeff*CopyNumber,
           CopyNumberFit=(GcRatioMedian-cnFitIntercept)/cnFitCoeff)

cnGcRatioPlot = ggplot(cnGcRatioSegments) +
  geom_point(aes(x=CopyNumber,y=CopyNumberFit),color='blue') +
  geom_abline(slope = 1,intercept = 0,color='red') +
  scale_x_log10() + scale_y_log10() +
  labs(x='Copy Number',y='GC Ratio Implied Copy Number')

ggsave(filename = paste0(outputDir, "/", sampleId, ".cn_gc_ratio_fit.png"), cnGcRatioPlot, units="cm",height=10,width=10,scale=1)
