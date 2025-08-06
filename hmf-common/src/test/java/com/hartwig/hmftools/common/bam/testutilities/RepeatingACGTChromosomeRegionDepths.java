package com.hartwig.hmftools.common.bam.testutilities;

import htsjdk.samtools.SAMFileWriter;

public class RepeatingACGTChromosomeRegionDepths extends ChromosomeRegionDepths
{
    public RepeatingACGTChromosomeRegionDepths(final int mChromosome)
    {
        super(mChromosome);
    }

    @Override
    public void writeEntriesForRange(final RegionDepth range, final SAMFileWriter bamWriter)
    {
        range.length100ReadsBamRegionWriter().writeEntries(bamWriter, ChromosomeWindow::toACGTBasesRegion);
    }
}
