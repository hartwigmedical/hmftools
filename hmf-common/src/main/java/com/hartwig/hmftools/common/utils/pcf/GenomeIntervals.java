package com.hartwig.hmftools.common.utils.pcf;

import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class GenomeIntervals
{
    public final List<PCFIntervals> mIntervals;

    public GenomeIntervals(final List<PCFIntervals> intervals)
    {
        PCFIntervals previousInterval = null;
        for(final PCFIntervals interval : intervals)
        {
            if(previousInterval != null)
            {
                Preconditions.checkArgument(previousInterval.mChromosome().compareTo(interval.mChromosome()) < 0);
            }
            previousInterval = interval;
        }
        this.mIntervals = intervals;
    }

    public List<ChrBaseRegion> regionsList()
    {
        List<ChrBaseRegion> regions = new ArrayList<>();
        mIntervals.forEach(intervals -> regions.addAll(intervals.mIntervals()));
        return regions;
    }
}
