package com.hartwig.hmftools.purple.region;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.google.common.collect.TreeMultimap;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.purple.segment.PurpleSupportSegment;
import com.hartwig.hmftools.purple.segment.SegmentRefiner;

public class ExcludedRegionsRefiner
{
    private final Map<String, SegmentRefiner> mRefiners = new HashMap<>();

    public ExcludedRegionsRefiner(Collection<ChrBaseRegion> regions)
    {
        TreeMultimap<String, GenomeRegion> byChr = TreeMultimap.create();
        regions.forEach(region -> byChr.put(region.chromosome(), region.genomeRegion()));

        byChr.asMap().forEach((chr, regionsOnChr) -> mRefiners.put(chr, new SegmentRefiner(new ArrayList<>(regionsOnChr))));
    }

    public boolean hasExcludedRegions(String chromosome)
    {
        return mRefiners.containsKey(chromosome);
    }

    public SegmentRefiner refiner(String chromosome)
    {
        return mRefiners.get(chromosome);
    }

    public List<PurpleSupportSegment> refine(List<PurpleSupportSegment> segments)
    {
        List<PurpleSupportSegment> result = new ArrayList<>();
        String currentChr = null;
        List<PurpleSupportSegment> currentSegments = null;
        for(PurpleSupportSegment segment : segments)
        {
            if(!segment.chromosome().equals(currentChr))
            {
                // Refine the segments and switch to the new chromosome
                result.addAll(refineList(currentChr, currentSegments));
                currentSegments = new ArrayList<>();
                currentChr = segment.chromosome();
            }
            currentSegments.add(segment);
        }
        result.addAll(refineList(currentChr, currentSegments));
        return result;
    }

    private List<PurpleSupportSegment> refineList(String chromosome, List<PurpleSupportSegment> segments)
    {
        if(chromosome == null)
        {
            return new ArrayList<>();
        }
        if(!hasExcludedRegions(chromosome))
        {
            return segments;
        }
        return refiner(chromosome).refine(segments);
    }
}
