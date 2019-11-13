package com.hartwig.hmftools.sage.context;

import java.util.Arrays;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class RefSequence {

    private static final int BUFFER = 1000;

    private int actualStart;
    private final ReferenceSequence sequence;

    public RefSequence(@NotNull final GenomeRegion region, @NotNull final IndexedFastaSequenceFile refGenome) {
        actualStart = Math.max(0, (int) region.start() - BUFFER);
        int maxLength = refGenome.getSequenceDictionary().getSequence(region.chromosome()).getSequenceLength();
        final int actualEnd = Math.min(maxLength, (int) region.end() + BUFFER);
        this.sequence = refGenome.getSubsequenceAt(region.chromosome(), actualStart, actualEnd);
    }

    public byte[] alignment(int start, int end) {
        return Arrays.copyOfRange(sequence.getBases(), start - actualStart, end - actualStart + 1);
    }
}
