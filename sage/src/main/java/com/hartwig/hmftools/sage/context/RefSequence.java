package com.hartwig.hmftools.sage.context;

import java.util.Arrays;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.read.IndexedBases;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class RefSequence {

    private static final int BUFFER = 1000;
    private final int sequenceLength;
    private final int actualStart;
    private final ReferenceSequence sequence;

    public RefSequence(@NotNull final GenomeRegion region, @NotNull final IndexedFastaSequenceFile refGenome) {
        sequenceLength = refGenome.getSequenceDictionary().getSequence(region.chromosome()).getSequenceLength();
        actualStart = Math.max(1, (int) region.start() - BUFFER);
        final int actualEnd = Math.min(sequenceLength, (int) region.end() + BUFFER);
        this.sequence = refGenome.getSubsequenceAt(region.chromosome(), actualStart, actualEnd);
    }

    private byte[] alignment(int start, int end) {
        return Arrays.copyOfRange(sequence.getBases(), start - actualStart, end - actualStart + 1);
    }

    /**
     * returns reference sequence spanning read with index pointing to alignment start
     */
    @NotNull
    public IndexedBases alignment(final SAMRecord record) {

        int alignmentStart = record.getAlignmentStart();

        int refPositionStart = record.getCigar().isLeftClipped()
                ? Math.max(1, alignmentStart - record.getCigar().getCigarElement(0).getLength())
                : alignmentStart;

        int alignmentEnd = Math.max(record.getAlignmentEnd(), alignmentStart + record.getReadLength());
        int alignmentStartIndex = alignmentStart - refPositionStart;

        return new IndexedBases(alignmentStart, alignmentStartIndex, alignment(refPositionStart, Math.min(alignmentEnd, sequenceLength)));
    }

}
