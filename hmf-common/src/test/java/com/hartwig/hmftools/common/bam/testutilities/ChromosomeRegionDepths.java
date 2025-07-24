package com.hartwig.hmftools.common.bam.testutilities;

import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Preconditions;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import htsjdk.samtools.SAMFileWriter;

public class ChromosomeRegionDepths
{
    public final int mChromosome;
    private final List<RegionDepth> ranges = new ArrayList<>();

    public ChromosomeRegionDepths(final int mChromosome)
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

    public void writeToBam(SAMFileWriter bamWriter, RefGenomeSource refGenomeSource)
    {
        ranges.forEach(range -> range.length100ReadsBamRegionWriter().writeEntries(bamWriter, refGenomeSource));
    }
}
