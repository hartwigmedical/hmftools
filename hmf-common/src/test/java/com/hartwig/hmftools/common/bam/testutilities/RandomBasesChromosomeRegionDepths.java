package com.hartwig.hmftools.common.bam.testutilities;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import htsjdk.samtools.SAMFileWriter;

public class RandomBasesChromosomeRegionDepths extends ChromosomeRegionDepths
{
    public RandomBasesChromosomeRegionDepths(final int mChromosome)
    {
        super(mChromosome);
    }

    @Override
    public void writeEntriesForRange(final RegionDepth range, final SAMFileWriter bamWriter)
    {
        range.length100ReadsBamRegionWriter().writeEntries(bamWriter, ChromosomeWindow::toRandomBasesRegion);
    }
}
