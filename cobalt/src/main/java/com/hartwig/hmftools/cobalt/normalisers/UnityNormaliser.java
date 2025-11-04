package com.hartwig.hmftools.cobalt.normalisers;

import com.hartwig.hmftools.cobalt.calculations.BamRatio;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class UnityNormaliser implements ResultsNormaliser
{
    private final DescriptiveStatistics RatiosStatistics = new DescriptiveStatistics();
    private final DescriptiveStatistics MegaBaseScaleRatioStatistics = new DescriptiveStatistics();
    private double Mean;
    private double MegaBaseScaleMean;

    @Override
    public void recordValue(final BamRatio bamRatio)
    {
        if(!bamRatio.mChromosome.isAutosome())
        {
            return;
        }
        double ratio = bamRatio.ratio();
        if(ratio > 0.0)
        {
            RatiosStatistics.addValue(ratio);
        }
        double mbScaleRatio = bamRatio.getDiploidAdjustedRatio();
        if(mbScaleRatio > 0.0)
        {
            MegaBaseScaleRatioStatistics.addValue(mbScaleRatio);
        }
    }

    @Override
    public void dataCollectionFinished()
    {
        Mean = RatiosStatistics.getMean();
        MegaBaseScaleMean = MegaBaseScaleRatioStatistics.getMean();
    }

    @Override
    public void normalise(BamRatio bamRatio)
    {
        bamRatio.normaliseByMean(Mean);
        if(MegaBaseScaleRatioStatistics.getN() > 0)
        {
            bamRatio.normaliseDiploidAdjustedRatio(MegaBaseScaleMean);
        }
    }
}
