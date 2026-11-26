package com.hartwig.hmftools.tars.liftback;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Genomic contamination zones (rRNA / 7SL / acrocentric) loaded from a TSV of Chromosome/PosStart/PosEnd.
public final class ExcludedRegions
{
    private final Map<String, List<BaseRegion>> mRegionsByChromosome;

    ExcludedRegions(final Map<String, List<BaseRegion>> regionsByChromosome)
    {
        mRegionsByChromosome = new HashMap<>();
        for(Map.Entry<String, List<BaseRegion>> entry : regionsByChromosome.entrySet())
        {
            List<BaseRegion> regions = new ArrayList<>(entry.getValue());
            Collections.sort(regions);
            mRegionsByChromosome.put(entry.getKey(), regions);
        }
    }

    public static ExcludedRegions load(final String filename)
    {
        Map<String, List<BaseRegion>> regions = ChrBaseRegion.loadChrBaseRegions(filename, false);
        if(regions == null)
        {
            throw new IllegalStateException(String.format("failed to read excluded regions file %s", filename));
        }
        return new ExcludedRegions(regions);
    }

    public boolean excludes(final String chromosome, final int posStart, final int posEnd)
    {
        List<BaseRegion> regions = mRegionsByChromosome.get(chromosome);
        if(regions == null)
        {
            return false;
        }

        // last region starting at or before posEnd; the start<=posEnd guard covers the "all regions after" case.
        BaseRegion region = regions.get(BaseRegion.binarySearch(posEnd, regions));
        return region.start() <= posEnd && region.end() >= posStart;
    }
}
