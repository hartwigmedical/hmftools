package com.hartwig.hmftools.cobalt.calculations;

import static com.hartwig.hmftools.cobalt.CobaltConstants.DEFAULT_GC_RATIO_MAX;
import static com.hartwig.hmftools.cobalt.CobaltConstants.DEFAULT_GC_RATIO_MIN;

import com.google.common.base.Preconditions;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class ReadDepthStatisticsNormaliser implements ResultsNormaliser
{
    private final DescriptiveStatistics descriptiveStatistics = new DescriptiveStatistics();
    private double Ratio = -1.0;
    private double mReadDepthMedian = -1.0;
    private double mReadDepthMean = -1.0;

    @Override
    public void dataCollectionFinished()
    {
        mReadDepthMean = descriptiveStatistics.getMean();
        mReadDepthMedian = descriptiveStatistics.getPercentile(50);
        Ratio =  mReadDepthMean / mReadDepthMedian;
    }

    public double readDepthMean()
    {
        Preconditions.checkState(Ratio > 0.0);
        return mReadDepthMean;
    }

    public double readDepthMedian()
    {
        Preconditions.checkState(Ratio > 0.0);
        return mReadDepthMedian;
    }

    @Override
    public void recordValue(final BamRatio bamRatio)
    {
        Preconditions.checkState(Ratio < 0.0);

        if (!bamRatio.mChromosome.isAutosome())
        {
            return;
        }
        if (bamRatio.gcContent() < DEFAULT_GC_RATIO_MIN)
        {
            return;
        }
        if(bamRatio.gcContent() > DEFAULT_GC_RATIO_MAX)
        {
            return;
        }
        double ratio = bamRatio.ratio();
        if (ratio > 0.0)
        {
            descriptiveStatistics.addValue(bamRatio.readDepth());
        }
    }

    @Override
    public void normalise(BamRatio bamRatio)
    {
        Preconditions.checkState(Ratio > 0.0);
        bamRatio.normaliseByMean(Ratio);
    }
}
