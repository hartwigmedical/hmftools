package com.hartwig.hmftools.common.utils.pcf;

import static com.google.common.base.Preconditions.checkArgument;

import java.util.List;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

public record PCFIntervals(HumanChromosome mChromosome, List<ChrBaseRegion> mIntervals)
{
    public PCFIntervals
    {
        ChrBaseRegion previousRegion = null;
        for(ChrBaseRegion region : mIntervals)
        {
            checkArgument(region.humanChromosome().equals(mChromosome));
            if(previousRegion != null)
            {
                checkArgument(region.start() > previousRegion.end());
            }
            previousRegion = region;
        }
    }
}
