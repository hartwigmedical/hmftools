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
somaticFile = args[4]
outputDir <- args[5]

print(sprintf('Producing Somatic plot for patient(%s) sample(%s), writing to output(%s)', patientId, sampleId, outputDir))

if (!file.exists(summaryFile))
{
  print(sprintf('Missing sample summary file: %s', summaryFile))
  stop()
}

if (!file.exists(somaticFile))
{
  print(sprintf('Missing sample somatic variants file: %s', cnSegmentFile))
  stop()
}

sampleSummary = read.csv(summaryFile,sep='\t') %>% filter(SampleId==sampleId)
allVariants = read.csv(somaticFile,sep='\t') %>% filter(SampleId==sampleId)

if(nrow(sampleSummary) == 0 | nrow(allVariants) == 0)
{
  print(sprintf('Empty filtered input files'))
  stop()
}

filteredVars = allVariants %>% filter(Filter=='PASS')

filteredVars = filteredVars %>% mutate(SampleVarVaf=ifelse(SampleDP>0,round(SampleAD/SampleDP,6),0))

# filtering and plotting threshold
qualPerAdThreshold=18
vafBucket=0.002
minVariantCount=5
minSampleDepth=20
minSampleDepthPerc=0.1

estimatedVaf=sampleSummary$SampleVAF
rawSomaticPurity=sampleSummary$RawSomaticPurity
peakPurity=sampleSummary$SNVPurity
weightedAvgDepth=sampleSummary$WeightedAvgDepth

depthThreshold=max(minSampleDepth,minSampleDepthPerc*weightedAvgDepth)

peakVars = filteredVars %>% filter(SampleDP>=depthThreshold&SampleQualPerAD>qualPerAdThreshold&SampleAD>0)

sampleAvgVaf = mean(peakVars$SampleVarVaf)

densityBandwidth = pmax(sampleAvgVaf/8, min(sampleAvgVaf/2, 0.01))

somaticPlot = ggplot() +
  geom_bar(data=peakVars %>% group_by(VafBucket=vafBucket*round(SampleVarVaf/vafBucket)) %>% count,stat="identity",position='identity',aes(y=n,x=VafBucket)) +
  geom_density(data=peakVars,aes(SampleVarVaf)) +
  geom_vline(xintercept=rawSomaticPurity*0.5,color='red') +
  geom_hline(yintercept=minVariantCount,color='blue') +
  labs(x='Filter Variant VAF',y='# Variants',title=sprintf('%s %s', patientId, sampleId))

if(peakPurity > rawSomaticPurity)
{
    somaticPlot = somaticPlot + geom_vline(xintercept=peakPurity*0.5,color='green')
}

ggsave(filename = paste0(outputDir, "/", sampleId, ".somatic_vaf.png"), somaticPlot, units="cm",height=10,width=10,scale=1)
