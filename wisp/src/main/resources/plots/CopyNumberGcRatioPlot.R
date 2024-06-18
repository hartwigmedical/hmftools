library(dplyr)
library(tidyr)
library(ggplot2)
library(gtable)

# Parse and check inputs
args <- commandArgs(trailingOnly = TRUE)
print(args)

if(length(args) < 5)
{
  print("Requires arguments 1=PatientId 2=SampleId 3=Plot calcs file 4=CN Segment file 5=Sample plot output directory")
  stop()
}

patientId = args[1]
sampleId = args[2]
plotCalcFile = args[3]
cnSegmentFile = args[4]
outputDir <- args[5]

print(sprintf('Producing CN-GC Ratio plot for patient(%s) sample(%s), writing to output(%s)', patientId, sampleId, outputDir))

if(!file.exists(plotCalcFile))
{
  print(sprintf('Missing sample summary file: %s', plotCalcFile))
  stop()
}

if(!file.exists(cnSegmentFile))
{
  print(sprintf('Missing sample CN-GcRatio file: %s', cnSegmentFile))
  stop()
}

samplePlotCalcs = read.csv(plotCalcFile,sep='\t')
cnData = read.csv(cnSegmentFile,sep='\t')

if('SampleId' %in% colnames(samplePlotCalcs))
{
    samplePlotCalcs = samplePlotCalcs %>% filter(SampleId==sampleId)
    cnData = cnData %>% filter(SampleId==sampleId)
}

if('PatientId' %in% colnames(samplePlotCalcs))
{
    samplePlotCalcs = samplePlotCalcs %>% filter(PatientId==patientId)
    cnData = cnData %>% filter(PatientId==patientId)
}

if(nrow(samplePlotCalcs) == 1 & nrow(cnData) > 0)
{
    purityEstimate=samplePlotCalcs$EstimatedPurity
    fitIntercept=samplePlotCalcs$FitIntercept
    fitCoefficient=samplePlotCalcs$FitCoefficient

    fitInterceptLow=samplePlotCalcs$FitInterceptLow
    fitCoefficientLow=samplePlotCalcs$FitCoefficientLow
    fitInterceptHigh=samplePlotCalcs$FitInterceptHigh
    fitCoefficientHigh=samplePlotCalcs$FitCoefficientHigh

    # thresholds for display
    cnLimitMax = 6
    cnLimitMin = 0
    cnLevelMargin = 0.2
    cnDataSizeFactor=0.003
    cnLevelSizeFactor=0.001
    gcRatioLimit = 2

    cnData = cnData %>% mutate(CnLevel=round(CopyNumber),
                               CnLevelDiff=abs(CopyNumber-CnLevel),
                               CnValid=CnLevel>=0&CnLevel<=6&CnLevelDiff<=cnLevelMargin,
                               SegmentWeight=sqrt(GcRatioCount),
                               GcRatioCapped=pmin(GcRatioMean,gcRatioLimit))

    cnLevelData = cnData %>% filter(CnValid) %>% group_by(CnLevel) %>%
      summarise(GcRatioWeightTotal=sum(SegmentWeight),
                GcRatioWeighted=pmin(sum(GcRatioMean*SegmentWeight)/GcRatioWeightTotal,gcRatioLimit))

    format_purity<-function(purity)
    {
      if(purity >= 0.01)
      {
        return (sprintf('%.4f', purity))
      } else {
        return (sprintf('%.6f', purity))
      }
    }

    plotTitle=sprintf('%s - %s: TF estimate(%s)', patientId, sampleId, format_purity(purityEstimate))

    titleSize=8
    labelSize=6
    alphaFade=0.75
    colourCnLevel='orange'
    colourCnRaw='blue'

    cnGcRatioPlot = ggplot() +
            geom_point(data=cnData %>% filter(CnLevel<=cnLimitMax),aes(x=CopyNumber,y=GcRatioCapped,size=SegmentWeight*cnDataSizeFactor),color=colourCnLevel) +
            geom_point(data=cnLevelData,aes(x=CnLevel,y=GcRatioWeighted,size=GcRatioWeightTotal*cnLevelSizeFactor),color=colourCnRaw,alpha=alphaFade) +
            scale_size_identity() +
            geom_abline(slope=fitCoefficient,intercept=fitIntercept,color='red') +
            geom_abline(slope=fitCoefficientLow,intercept=fitInterceptLow,color='grey') +
            geom_abline(slope=fitCoefficientHigh,intercept=fitInterceptHigh,color='grey') +
            theme(legend.position="none") +
            labs(x='Tumor Copy Number',y='Sample GC Ratio',title=plotTitle) +
            theme(plot.title=element_text(size=titleSize),
                  axis.title=element_text(size=labelSize),axis.text=element_text(size=labelSize),
                  legend.position="none")

    ggsave(filename = paste0(outputDir, "/", sampleId, ".cn_gc_ratio_fit.png"), cnGcRatioPlot, units="cm",height=10,width=10,scale=1)

} else
{
  print(sprintf('Empty filtered input files'))
}
