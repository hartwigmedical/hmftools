package com.hartwig.hmftools.cobalt.calculations;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class BamRatiosMeanCalculator
{
    private final DescriptiveStatistics descriptiveStatistics = new DescriptiveStatistics();

    public void recordValue(BamRatio bamRatio)
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

    public double mean()
    {
        return descriptiveStatistics.getMean();
    }
}
