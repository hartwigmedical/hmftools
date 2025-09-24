package com.hartwig.hmftools.cobalt.calculations;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class UnityNormaliser implements ResultsNormaliser
{
    private final DescriptiveStatistics descriptiveStatistics = new DescriptiveStatistics();

    @Override
    public void recordValue(final BamRatio bamRatio)
    {
        if (!bamRatio.mChromosome.isAutosome())
        {
            return;
        }
        double ratio = bamRatio.ratio();
        if (ratio > 0.0)
        {
            descriptiveStatistics.addValue(ratio);
        }
    }

    @Override
    public void applyNormalisation(BamRatio bamRatio)
    {
        double mean = descriptiveStatistics.getMean();
        bamRatio.normaliseByMean(mean);
    }
}
