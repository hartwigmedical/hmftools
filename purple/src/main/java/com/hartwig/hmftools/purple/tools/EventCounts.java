package com.hartwig.hmftools.purple.tools;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.SortedMap;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

public class EventCounts
{
    private final SortedMap<HumanChromosome, ChromosomeRegionCounts> mChromosomeRegionCountsMap;

    public EventCounts(final SortedMap<HumanChromosome, ChromosomeRegionCounts> mChromosomeRegionCountsMap)
    {
        this.mChromosomeRegionCountsMap = mChromosomeRegionCountsMap;
    }

    public void mergeAdd(EventCounts other)
    {
        for(HumanChromosome chromosome : other.mChromosomeRegionCountsMap.keySet())
        {
            if(mChromosomeRegionCountsMap.containsKey(chromosome))
            {
                ChromosomeRegionCounts existing = mChromosomeRegionCountsMap.get(chromosome);
                existing.addAll(other.mChromosomeRegionCountsMap.get(chromosome));
            }
            else
            {
                mChromosomeRegionCountsMap.put(chromosome, other.mChromosomeRegionCountsMap.get(chromosome));
            }
        }
    }

    public SortedMap<HumanChromosome, ChromosomeRegionCounts> data()
    {
        return mChromosomeRegionCountsMap;
    }

    public void writeTo(BufferedWriter bw, final int minCountForInclusion, final RefGenomeVersion refGenomeVersion) throws IOException
    {
        bw.write("Chromosome,RegionStart,RegionEnd,Type,Frequency\n");
        for(HumanChromosome chromosome : mChromosomeRegionCountsMap.keySet())
        {
            String chr = refGenomeVersion.versionedChromosome(chromosome);
            SortedMap<RegionAmpDel, Integer> countMap = mChromosomeRegionCountsMap.get(chromosome).counts();
            for(RegionAmpDel re : countMap.keySet())
            {
                int count = countMap.get(re);
                if(count >= minCountForInclusion)
                {
                    String line = String.format("%s,%d,%d,%s,%d\n", chr, re.start(), re.end(), re.type(), count);
                    bw.write(line);
                }
            }
        }
        bw.flush();
    }
}

