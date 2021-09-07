package com.hartwig.hmftools.sage.ref;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;
import com.hartwig.hmftools.sage.read.IndexedBases;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class RefSequence
{
    public final int End;
    public final int Start;
    public final ReferenceSequence Sequence;
    public final IndexedBases IndexedBases;

    private static final int BUFFER = 1000;

    public RefSequence(final ChrBaseRegion region, final ReferenceSequenceFile refGenome)
    {
        final int sequenceEnd = refGenome.getSequenceDictionary().getSequence(region.Chromosome).getSequenceLength();
        Start = Math.max(1, region.start() - BUFFER);
        End = Math.min(sequenceEnd, region.end() + BUFFER);
        Sequence = refGenome.getSubsequenceAt(region.Chromosome, Start, End);
        IndexedBases = new IndexedBases(Start, 0, Sequence.getBases());
    }

    @VisibleForTesting
    public RefSequence(final ReferenceSequence sequence)
    {
        Sequence = sequence;
        Start = sequence.getContigIndex() + 1;
        End = Start + sequence.getBases().length - 1;
        IndexedBases = new IndexedBases(Start, 0, sequence.getBases());
    }

    // returns reference sequence spanning read with index pointing to alignment start
    public IndexedBases alignment()
    {
        return IndexedBases;
    }

}
