package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.tars.common.TarsConstants.TARS_LOGGER;

import java.util.Comparator;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

// Genomic contamination zones (rRNA / 7SL / acrocentric) loaded from a TSV (Chromosome/PosStart/PosEnd), queried
// on genomic coordinates. Per-chromosome, sorted for O(log n) overlap tests.
public final class ExcludedRegions
{
    private final Map<String, List<BaseRegion>> mRegionsByChromosome;

    ExcludedRegions(final Map<String, List<BaseRegion>> regionsByChromosome)
    {
        mRegionsByChromosome = regionsByChromosome;
        for(List<BaseRegion> regions : mRegionsByChromosome.values())
            regions.sort(Comparator.comparingInt(BaseRegion::start));
    }

    public static ExcludedRegions load(final String filename)
    {
        Map<String, List<BaseRegion>> regions = ChrBaseRegion.loadChrBaseRegions(filename, false);
        if(regions == null)
        {
            TARS_LOGGER.error("failed to read excluded regions file {}", filename);
            System.exit(1);
        }
        return new ExcludedRegions(regions);
    }

    public boolean excludes(final String chromosome, final int posStart, final int posEnd)
    {
        List<BaseRegion> regions = mRegionsByChromosome.get(chromosome);
        if(regions == null)
            return false;

        // last region starting at or before posEnd; the start<=posEnd guard covers the "all regions after" case.
        BaseRegion region = regions.get(BaseRegion.binarySearch(posEnd, regions));
        return region.start() <= posEnd && region.end() >= posStart;
    }
}
