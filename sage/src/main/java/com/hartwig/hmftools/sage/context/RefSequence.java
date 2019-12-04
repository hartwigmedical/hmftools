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
        actualStart = Math.max(0, (int) region.start() - BUFFER);
        final int actualEnd = Math.min(sequenceLength, (int) region.end() + BUFFER);
        this.sequence = refGenome.getSubsequenceAt(region.chromosome(), actualStart, actualEnd);
    }

    private byte[] alignment(int start, int end) {
        return Arrays.copyOfRange(sequence.getBases(), start - actualStart, end - actualStart + 1);
    }

    @NotNull
    public IndexedBases alignment(final SAMRecord record) {

        int alignmentStart = record.getAlignmentStart();

        int readStart = alignmentStart;
//                record.getCigar().isLeftClipped() ? alignmentStart - record.getCigar().getCigarElement(0).getLength() : alignmentStart;

        int alignmentEnd = Math.max(record.getAlignmentEnd(), alignmentStart + record.getReadLength());
        int index = alignmentStart - readStart;

        return new IndexedBases(index, alignment(Math.max(readStart, 0), Math.min(alignmentEnd, sequenceLength)));
    }

}
