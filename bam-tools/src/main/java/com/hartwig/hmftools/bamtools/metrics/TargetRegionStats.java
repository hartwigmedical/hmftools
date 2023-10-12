package com.hartwig.hmftools.bamtools.metrics;

import static java.lang.String.format;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class TargetRegionStats
{
    public final ChrBaseRegion Region;
    public ReadCounts Counts;

    public TargetRegionStats(final ChrBaseRegion region)
    {
        Region = region;
        Counts = new ReadCounts();
    }

    public void mergeCounts(final TargetRegionStats other)
    {
        Counts.merge(other.Counts);
    }

    public String toString()
    {
        return format("%s reads(%s)", Region, Counts);
    }
}
