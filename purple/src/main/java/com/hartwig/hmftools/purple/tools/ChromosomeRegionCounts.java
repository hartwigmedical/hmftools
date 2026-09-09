package com.hartwig.hmftools.purple.tools;

import java.util.List;
import java.util.NavigableMap;
import java.util.SortedMap;
import java.util.TreeMap;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

public class ChromosomeRegionCounts
{
    private final Chromosome mChromosome;
    private final NavigableMap<RegionAmpDel, Integer> mEventCounts = new TreeMap<>();

    public ChromosomeRegionCounts(Chromosome chromosome, List<RegionGeneEvents> events)
    {
        mChromosome = chromosome;
        for(RegionGeneEvents event : events)
        {
            Preconditions.checkArgument(event.chromosome().equals(mChromosome));
            final RegionAmpDel ampDel = new RegionAmpDel(event);
            mEventCounts.merge(ampDel, 1, Integer::sum);
        }
    }

    public SortedMap<RegionAmpDel, Integer> counts()
    {
        return mEventCounts;
    }

    public void addAll(ChromosomeRegionCounts other)
    {
        Preconditions.checkArgument(other.mChromosome.equals(mChromosome));
        for(RegionAmpDel regionCount : other.mEventCounts.keySet())
        {
            mEventCounts.merge(regionCount, other.mEventCounts.get(regionCount), Integer::sum);
        }
    }
}
