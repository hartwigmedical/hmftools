package com.hartwig.hmftools.sage.common;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class RefSequence
{
    public final int End;
    public final int Start;
    public final IndexedBases IndexedBases;

    private static final int BUFFER = 1000;

    public RefSequence(final ChrBaseRegion region, final RefGenomeInterface refGenome)
    {
        final int sequenceEnd = refGenome.getChromosomeLength(region.Chromosome);
        Start = Math.max(1, region.start() - BUFFER);
        End = Math.min(sequenceEnd, region.end() + BUFFER);
        final byte[] seqBases = refGenome.getBases(region.Chromosome, Start, End);
        IndexedBases = new IndexedBases(Start, 0, seqBases);
    }

    @VisibleForTesting
    public RefSequence(final ChrBaseRegion region, final ReferenceSequenceFile refGenome)
    {
        final int sequenceEnd = refGenome.getSequenceDictionary().getSequence(region.Chromosome).getSequenceLength();
        Start = Math.max(1, region.start() - BUFFER);
        End = Math.min(sequenceEnd, region.end() + BUFFER);
        final byte[] seqBases = refGenome.getSubsequenceAt(region.Chromosome, Start, End).getBases();
        IndexedBases = new IndexedBases(Start, 0, seqBases);
    }

    @VisibleForTesting
    public RefSequence(final ReferenceSequence sequence)
    {
        Start = sequence.getContigIndex() + 1;
        End = Start + sequence.getBases().length - 1;
        IndexedBases = new IndexedBases(Start, 0, sequence.getBases());
    }

    // returns reference sequence spanning read with index pointing to alignment start
    public IndexedBases alignment() { return IndexedBases; }

}
