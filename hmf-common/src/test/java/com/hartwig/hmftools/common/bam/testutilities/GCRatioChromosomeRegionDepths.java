package com.hartwig.hmftools.common.bam.testutilities;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import htsjdk.samtools.SAMFileWriter;

public class GCRatioChromosomeRegionDepths extends ChromosomeRegionDepths
{
    private final int gcPercentage;

    public GCRatioChromosomeRegionDepths(final HumanChromosome mChromosome, final int gcPercentage)
    {
        super(mChromosome);
        this.gcPercentage = gcPercentage;
    }

    @Override
    public void writeEntriesForRange(final RegionDepth range, final SAMFileWriter bamWriter)
    {
        range.length100ReadsBamRegionWriter().writeEntries(bamWriter, window -> window.toBasesWithGivenGC(gcPercentage));
    }
}
