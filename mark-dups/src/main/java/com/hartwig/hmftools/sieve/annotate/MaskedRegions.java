package com.hartwig.hmftools.sieve.annotate;

import java.util.HashMap;
import java.util.List;
import java.util.TreeSet;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class MaskedRegions
{
    private final HashMap<String, TreeSet<ChrBaseRegion>> mSortedRegionsByChr;

    public MaskedRegions(final List<ChrBaseRegion> regions)
    {
        mSortedRegionsByChr = new HashMap<>();
        for(ChrBaseRegion region : regions)
        {
            if(!mSortedRegionsByChr.containsKey(region.Chromosome))
            {
                mSortedRegionsByChr.put(region.Chromosome, new TreeSet<>(
                        (final ChrBaseRegion q1, final ChrBaseRegion q2) ->
                        {
                            final int startDiff = q1.start() - q2.start();
                            if(startDiff != 0)
                            {
                                return startDiff;
                            }

                            return q1.end() - q2.end();
                        }));
            }

            mSortedRegionsByChr.get(region.Chromosome).add(region);
        }
    }

    int distance(final ChrBaseRegion region)
    {
        if(!mSortedRegionsByChr.containsKey(region.Chromosome))
        {
            return -1;
        }

        final TreeSet<ChrBaseRegion> sortedRegions = mSortedRegionsByChr.get(region.Chromosome);
        final ChrBaseRegion floor = sortedRegions.floor(region);
        final ChrBaseRegion ceiling = sortedRegions.ceiling(region);

        int floorDist = -1;
        int ceilingDist = -1;
        if(floor != null)
        {
            floorDist = (floor.end() <= region.start()) ? region.start() - floor.end() : 0;
        }

        if(ceiling != null)
        {
            ceilingDist = (region.end() <= ceiling.start()) ? ceiling.start() - region.end() : 0;
        }

        if(floorDist == -1)
        {
            return ceilingDist;
        }

        if(ceilingDist == -1)
        {
            return floorDist;
        }

        return Math.min(floorDist, ceilingDist);
    }
}
