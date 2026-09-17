package com.hartwig.hmftools.purple.tools;

import java.util.Comparator;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class EventRegionComparator implements Comparator<ChrBaseRegion>
{
    @Override
    public int compare(final ChrBaseRegion o1, final ChrBaseRegion o2)
    {
        int result = o1.humanChromosome().compareTo(o2.humanChromosome());
        if(result != 0)
        {
            return result;
        }
        result = o1.start() - o2.start();
        if(result != 0)
        {
            return result;
        }
        result = o1.end() - o2.end();
        return result;
    }
}
