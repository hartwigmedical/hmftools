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
variantImpliedTFs = read.csv(somaticPeakFile,sep='\t')

if('SampleId' %in% colnames(sampleSummary))
{
    sampleSummary = sampleSummary %>% filter(SampleId==sampleId)
    variantImpliedTFs = variantImpliedTFs %>% filter(SampleId==sampleId)
}

if('PatientId' %in% colnames(sampleSummary))
{
    sampleSummary = sampleSummary %>% filter(PatientId==patientId)
    variantImpliedTFs = variantImpliedTFs %>% filter(PatientId==patientId)
}

clonalMethod = 'NONE'

if(nrow(sampleSummary) == 1)
{
    clonalMethod=sampleSummary$ClonalMethod
}

format_purity<-function(purity)
{
    if(purity >= 0.01)
    {
        return (sprintf('%.4f', purity))
    } else {
        return (sprintf('%.6f', purity))
    }
}

if(nrow(sampleSummary) == 1 & nrow(variantImpliedTFs) > 0 & clonalMethod != 'NONE')
{
    # filtering and plotting threshold
    minVariantCount = 3

    rawSomaticPurity = sampleSummary$RawSNVPurity
    maxImpliedTF = max(variantImpliedTFs$ImpliedTF)
    impliedTfBucket = pmin(0.01,maxImpliedTF/100)

    variantImpliedTFs = variantImpliedTFs %>% mutate(ImpliedTF=pmin(ImpliedTF,2.0))

    peakPurity = sampleSummary$SNVPurity
    peakPurityLow = sampleSummary$SNVPurityLow
    peakPurityHigh = sampleSummary$SNVPurityHigh
    densityBandwidth = sampleSummary$PeakBandwidth
    densityBandwidthLow = sampleSummary$PeakBandwidthLow
    densityBandwidthHigh = sampleSummary$PeakBandwidthHigh

    # colours used in plot
    colourKdeLine = 'orange2' # black
    colourFreqDist = 'grey'
    colourRawPurity = 'red'
    colourPeakPurity = 'orange2'
    colourLowPurity = 'grey'
    colourHighPurity = 'grey'
    titleSize=8
    labelSize=6

    if(peakPurity > rawSomaticPurity)
    {
        purityStr = sprintf('TF estimate(%s raw=%s)', format_purity(peakPurity), format_purity(rawSomaticPurity))
    } else {
        purityStr = sprintf('TF estimate(%s)', format_purity(rawSomaticPurity))
    }

    plotTitle=sprintf('%s - %s: %s, %s bandwidth(%.4f)',patientId,sampleId,purityStr,clonalMethod,densityBandwidth)

    somaticPlot = ggplot() +
      geom_bar(data=variantImpliedTFs %>% group_by(ImpliedTfBucket=impliedTfBucket*round(ImpliedTF/impliedTfBucket)) %>% count,
      stat="identity",position='identity',aes(y=n,x=ImpliedTfBucket),fill=colourFreqDist) +
      geom_vline(xintercept=rawSomaticPurity,color=colourRawPurity) +
      labs(x='Variant Implied TF',y='# Variants',title=plotTitle) +
      theme(plot.title=element_text(size=titleSize),
            axis.title=element_text(size=labelSize),axis.text=element_text(size=labelSize))

    if(clonalMethod == 'VAF_PEAK' | clonalMethod == 'NO_PEAK')
    {
        somaticPlot = somaticPlot +
            geom_density(data=variantImpliedTFs,bw=densityBandwidth,mapping = aes(x=ImpliedTF,y=after_stat(scaled)*20),color=colourKdeLine)

        if(densityBandwidthLow != densityBandwidth)
        {
            somaticPlot = somaticPlot +
                geom_density(data=variantImpliedTFs,bw=densityBandwidthLow,mapping = aes(x=ImpliedTF,y=after_stat(scaled)*20),color=colourLowPurity)
        }

        if(densityBandwidthHigh != densityBandwidth)
        {
            somaticPlot = somaticPlot +
                geom_density(data=variantImpliedTFs,bw=densityBandwidthHigh,mapping = aes(x=ImpliedTF,y=after_stat(scaled)*20),color=colourHighPurity)
        }
    }

    if(clonalMethod!='NO_PEAK')
    {
        peakTF = peakPurity
        somaticPlot = somaticPlot + geom_vline(xintercept=peakTF,color=colourPeakPurity)

        if(peakPurityLow > 0 & peakPurityLow < peakPurity)
        {
            peakTFLow = peakPurityLow
            somaticPlot = somaticPlot + geom_vline(xintercept=peakTFLow,color=colourLowPurity)
        }

        if(peakPurityHigh > 0 & peakPurityHigh > peakPurity)
        {
            peakTFHigh = peakPurityHigh
            somaticPlot = somaticPlot + geom_vline(xintercept=peakTFHigh,color=colourHighPurity)
        }
    }

    ggsave(filename = paste0(outputDir, "/", sampleId, ".somatic_vaf.png"), somaticPlot, units="cm",height=10,width=15,scale=1)

} else
{
  print(sprintf('Empty filtered input files'))
}
