package com.hartwig.hmftools.redux.splice;

import java.util.Comparator;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.redux.unmap.ReadUnmapper;

import htsjdk.samtools.SAMRecord;

// Genomic regions a fragment is dropped for before liftback (e.g. RNA rRNA / 7SL contamination zones). Loaded
// from a regions TSV (Chromosome/PosStart/PosEnd) and queried per fragment: if a mapped primary read overlaps
// an excluded region the whole fragment is filtered, so contaminating reads never get lifted or reach dedup.
// Per-chromosome regions are sorted for an O(log N) overlap test.
public class ExcludedRegions
{
    private final Map<String, List<ChrBaseRegion>> mRegionsByChr;

    ExcludedRegions(final Map<String, List<ChrBaseRegion>> regionsByChr)
    {
        mRegionsByChr = regionsByChr;
        for(final List<ChrBaseRegion> regions : mRegionsByChr.values())
            regions.sort(Comparator.comparingInt(ChrBaseRegion::start));
    }

    public static ExcludedRegions load(final String filename)
    {
        return new ExcludedRegions(ReadUnmapper.loadAlwaysUnmapRegions(filename));
    }

    // a fragment is excluded if any of its mapped primary reads overlaps an excluded region. Supplementary and
    // unmapped reads aren't tested -- the primary placement decides the fragment.
    public boolean fragmentExcluded(final List<SAMRecord> group)
    {
        for(final SAMRecord record : group)
        {
            if(record.getReadUnmappedFlag() || record.getSupplementaryAlignmentFlag() || record.isSecondaryAlignment())
                continue;

            if(overlaps(record.getReferenceName(), record.getAlignmentStart(), record.getAlignmentEnd()))
                return true;
        }
        return false;
    }

    private boolean overlaps(final String chromosome, final int posStart, final int posEnd)
    {
        final List<ChrBaseRegion> regions = mRegionsByChr.get(chromosome);
        if(regions == null || regions.isEmpty())
            return false;

        // rightmost region whose start <= posEnd; non-overlapping regions, so it's the only candidate.
        int lo = 0;
        int hi = regions.size() - 1;
        int candidate = -1;
        while(lo <= hi)
        {
            final int mid = (lo + hi) >>> 1;
            if(regions.get(mid).start() <= posEnd)
            {
                candidate = mid;
                lo = mid + 1;
            }
            else
            {
                hi = mid - 1;
            }
        }

        return candidate >= 0 && regions.get(candidate).end() >= posStart;
    }
}
