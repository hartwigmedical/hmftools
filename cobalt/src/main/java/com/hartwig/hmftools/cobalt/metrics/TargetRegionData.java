package com.hartwig.hmftools.cobalt.metrics;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class TargetRegionData extends ChrBaseRegion
{
    public FragmentGcMap FragmentGcCounts;
    private final DescriptiveStatistics Statistics = new DescriptiveStatistics();

//    public TargetRegionData(final String chromosome, final BaseRegion region)
//    {
//        super(chromosome, region.start(), region.end());
//
//        FragmentGcCounts = new FragmentGcMap();
//    }

    public TargetRegionData(final String chromosome, int start, int end)
    {
        super(chromosome, start, end);
    }

    public void recordFragment(int length)
    {
        Statistics.addValue(length);
    }

    public WindowStatistics statistics()
    {
        return new WindowStatistics(chromosome(), start(),
                Statistics.getN(), Statistics.getMean(),
                Statistics.getPercentile(50.0), Statistics.getMin(),
                Statistics.getMax(), Statistics.getStandardDeviation());
    }

    public String toString() { return super.toString(); }
}
