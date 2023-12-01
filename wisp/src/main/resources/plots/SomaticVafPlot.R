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

if (!file.exists(summaryFile))
{
  print(sprintf('Missing sample summary file: %s', summaryFile))
  stop()
}

if (!file.exists(somaticPeakFile))
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

if(nrow(sampleSummary) == 0 | nrow(variantVafRatios) == 0)
{
  print(sprintf('Empty filtered input files'))
  stop()
}

# no further filters to apppy
#filteredVars = allVariants %>% filter(Filter=='PASS')
#filteredVars = filteredVars %>% mutate(SampleVarVaf=ifelse(SampleDP>0,round(SampleAD/SampleDP,6),0))
#minSampleDepth=20
#minSampleDepthPerc=0.1

# filtering and plotting threshold
vafBucket=0.1
vafRatioMax=8
minVariantCount=5

adjSampleVaf=sampleSummary$AdjSampleVaf
rawSomaticPurity=sampleSummary$RawSomaticPurity
peakPurity=sampleSummary$SNVPurity
weightedAvgDepth=sampleSummary$WeightedAvgDepth
clonalMethod=sampleSummary$ClonalMethod

#peakVars = filteredVars %>% filter(SampleDP>=depthThreshold&SampleQualPerAD>qualPerAdThreshold&SampleAD>0)
#sampleAvgVaf = mean(peakVars$SampleVarVaf)
#densityBandwidth = pmax(sampleAvgVaf/8, min(sampleAvgVaf/2, 0.01))
densityBandwidth=0.15

rawVafRatio = 1

somaticPlot = ggplot() +
  geom_bar(data=variantVafRatios %>% group_by(VafRatioBucket=vafBucket*round(VafRatio/vafBucket)) %>% count,stat="identity",position='identity',aes(y=n,x=VafRatioBucket)) +
  geom_density(data=variantVafRatios,aes(VafRatio)) +
  geom_vline(xintercept=rawVafRatio,color='red') +
  geom_hline(yintercept=minVariantCount,color='blue') +
  labs(x='Variant VAF Ratio',y='# Variants',title=sprintf('Patient %s, Sample %s, %s',patientId,sampleId,clonalMethod))

if(clonalMethod!='NO_PEAK' & peakPurity > rawSomaticPurity)
{
    peakVafRatio = peakPurity/rawSomaticPurity
    somaticPlot = somaticPlot + geom_vline(xintercept=peakVafRatio  ,color='green')
}

ggsave(filename = paste0(outputDir, "/", sampleId, ".somatic_vaf.png"), somaticPlot, units="cm",height=10,width=10,scale=1)
