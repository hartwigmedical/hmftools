package com.hartwig.hmftools.sage.context;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.read.IndexedBases;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class RefSequence {

    private static final int BUFFER = 1000;
    private final int end;
    private final int start;
    private final ReferenceSequence sequence;

    public RefSequence(@NotNull final GenomeRegion region, @NotNull final IndexedFastaSequenceFile refGenome) {
        end = refGenome.getSequenceDictionary().getSequence(region.chromosome()).getSequenceLength();
        start = Math.max(1, (int) region.start() - BUFFER);
        final int actualEnd = Math.min(this.end, (int) region.end() + BUFFER);
        this.sequence = refGenome.getSubsequenceAt(region.chromosome(), start, actualEnd);
    }

    @VisibleForTesting
    RefSequence(final ReferenceSequence sequence) {
        this.sequence = sequence;
        this.start = sequence.getContigIndex() + 1;
        this.end = start + sequence.getBases().length - 1;
    }

    /**
     * returns reference sequence spanning read with index pointing to alignment start
     */
    @NotNull
    public IndexedBases alignment() {
        return new IndexedBases(start, 0, sequence.getBases());
    }

}
