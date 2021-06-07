package com.hartwig.hmftools.sage.ref;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.read.IndexedBases;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;

public class RefSequence
{

    private static final int BUFFER = 1000;

    private final int end;
    private final int start;
    private final ReferenceSequence sequence;
    private final IndexedBases indexedBases;

    public RefSequence(@NotNull final GenomeRegion region, @NotNull final ReferenceSequenceFile refGenome)
    {
        final int sequenceEnd = refGenome.getSequenceDictionary().getSequence(region.chromosome()).getSequenceLength();
        this.start = Math.max(1, (int) region.start() - BUFFER);
        this.end = Math.min(sequenceEnd, (int) region.end() + BUFFER);
        this.sequence = refGenome.getSubsequenceAt(region.chromosome(), start, end);
        this.indexedBases = new IndexedBases(start, 0, sequence.getBases());
    }

    @VisibleForTesting
    public RefSequence(final ReferenceSequence sequence)
    {
        this.sequence = sequence;
        this.start = sequence.getContigIndex() + 1;
        this.end = start + sequence.getBases().length - 1;
        this.indexedBases = new IndexedBases(start, 0, sequence.getBases());
    }

    /**
     * returns reference sequence spanning read with index pointing to alignment start
     */
    @NotNull
    public IndexedBases alignment()
    {
        return indexedBases;
    }

}
