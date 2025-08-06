package com.hartwig.hmftools.common.bam.testutilities;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import htsjdk.samtools.SAMFileWriter;

public class GenomeBackedChromosomeRegionDepths extends ChromosomeRegionDepths
{
    private final RefGenomeSource mRefGenome;

    public GenomeBackedChromosomeRegionDepths(final int mChromosome, final RefGenomeSource refGenomeSource)
    {
        super(mChromosome);
        this.mRefGenome = refGenomeSource;
    }

    @Override
    public void writeEntriesForRange(final RegionDepth range, final SAMFileWriter bamWriter)
    {
        range.length100ReadsBamRegionWriter().writeEntries(bamWriter, chromosomeWindow -> chromosomeWindow.toBaseRegionPair(mRefGenome));
    }
}
