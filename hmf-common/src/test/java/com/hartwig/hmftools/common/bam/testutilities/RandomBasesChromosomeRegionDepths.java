package com.hartwig.hmftools.common.bam.testutilities;

import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import htsjdk.samtools.SAMFileWriter;

public class RandomBasesChromosomeRegionDepths extends ChromosomeRegionDepths
{
    public RandomBasesChromosomeRegionDepths(final HumanChromosome mChromosome)
    {
        super(mChromosome);
    }

    @Override
    public void writeEntriesForRange(final RegionDepth range, final SAMFileWriter bamWriter)
    {
        range.length100ReadsBamRegionWriter().writeEntries(bamWriter, ChromosomeWindow::toRandomBasesRegion);
    }
}
