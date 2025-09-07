package com.hartwig.hmftools.common.bam.testutilities;

import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import htsjdk.samtools.SAMFileWriter;

public abstract class ChromosomeRegionDepths
{
    public final HumanChromosome mChromosome;
    private final List<RegionDepth> ranges = new ArrayList<>();

    public ChromosomeRegionDepths(final HumanChromosome mChromosome)
    {
        this.mChromosome = mChromosome;
    }

    public void addRange(final int start, final int end, final int depth)
    {
        RegionDepth newRange = new RegionDepth(mChromosome, start, end, depth);
        if(!ranges.isEmpty())
        {
            RegionDepth last = ranges.get(ranges.size() - 1);
            Preconditions.checkArgument(last.precedes(newRange));
        }
        ranges.add(newRange);
    }

    public void writeToBam(SAMFileWriter bamWriter)
    {
        ranges.forEach(range -> writeEntriesForRange(range, bamWriter));
    }

    public abstract void writeEntriesForRange(RegionDepth regionDepth, SAMFileWriter bamWriter);
}
