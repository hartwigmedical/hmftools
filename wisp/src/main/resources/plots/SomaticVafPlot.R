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
somaticPeakFile = args[4]
outputDir <- args[5]

print(sprintf('Producing Somatic plot for patient(%s) sample(%s), writing to output(%s)', patientId, sampleId, outputDir))

if(!file.exists(summaryFile))
{
  print(sprintf('Missing sample summary file: %s', summaryFile))
  stop()
}

if(!file.exists(somaticPeakFile))
{
  print(sprintf('Missing sample somatic variants file: %s', cnSegmentFile))
  stop()
}

sampleSummary = read.csv(summaryFile,sep='\t')
variantVafRatios = read.csv(somaticPeakFile,sep='\t')

if('SampleId' %in% colnames(sampleSummary))
{
    sampleSummary = sampleSummary %>% filter(SampleId==sampleId)
    variantVafRatios = variantVafRatios %>% filter(SampleId==sampleId)
}

clonalMethod=sampleSummary$ClonalMethod

if(nrow(sampleSummary) == 1 & nrow(variantVafRatios) > 0 & clonalMethod!='NONE')
{
    # no further filters to apppy

    # filtering and plotting threshold
    vafBucket=0.1
    vafRatioMax=8
    minVariantCount=5

    adjSampleVaf=sampleSummary$AdjSampleVaf
    rawSomaticPurity=sampleSummary$RawSomaticPurity
    peakPurity=sampleSummary$SNVPurity
    weightedAvgDepth=sampleSummary$WeightedAvgDepth
    tumorPurity=sampleSummary$TumorPurity

    ratioSummary = variantVafRatios %>%
      summarise(MeanRatio=mean(VafRatio),MedianRatio=median(VafRatio)) %>%
      mutate(MeanVsMedian=round(pmax(MeanRatio/MedianRatio,1),2))

    densityBandwidth=0.15*sqrt(ratioSummary$MeanVsMedian)

    rawVafRatio = 1

    plotTitle=sprintf('%s - %s: purity=%.2f, %s bandwidth(%.2f)',patientId,sampleId,tumorPurity,clonalMethod,densityBandwidth)

    somaticPlot = ggplot() +
      geom_bar(data=variantVafRatios %>% group_by(VafRatioBucket=vafBucket*round(VafRatio/vafBucket)) %>% count,
      stat="identity",position='identity',aes(y=n,x=VafRatioBucket),fill='grey') +
      geom_density(data=variantVafRatios,,mapping = aes(x=VafRatio,y=after_stat(scaled)*20),color='black') +
      geom_vline(xintercept=rawVafRatio,color='red') +
      geom_hline(yintercept=minVariantCount,color='blue') +
      labs(x='Variant VAF Ratio',y='# Variants',title=plotTitle) +
      theme(plot.title=element_text(size=8),
            axis.title=element_text(size=6),axis.text=element_text(size=6))

    if(clonalMethod!='NO_PEAK' & peakPurity > rawSomaticPurity)
    {
        peakVafRatio = peakPurity/rawSomaticPurity
        somaticPlot = somaticPlot + geom_vline(xintercept=peakVafRatio  ,color='green')
    }

    ggsave(filename = paste0(outputDir, "/", sampleId, ".somatic_vaf.png"), somaticPlot, units="cm",height=10,width=15,scale=1)

} else
{
  print(sprintf('Empty filtered input files'))
}
