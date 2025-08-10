package com.hartwig.hmftools.common.bam.testutilities;

import htsjdk.samtools.SAMFileWriter;

public class GCRatioChromosomeRegionDepths extends ChromosomeRegionDepths
{
    private final double gcRatio;

    public GCRatioChromosomeRegionDepths(final int mChromosome, final double gcRatio)
    {
        super(mChromosome);
        this.gcRatio = gcRatio;
    }

    @Override
    public void writeEntriesForRange(final RegionDepth range, final SAMFileWriter bamWriter)
    {
        range.length100ReadsBamRegionWriter().writeEntries(bamWriter, window -> window.toBasesWithGivenGC(gcRatio));
    }
}
