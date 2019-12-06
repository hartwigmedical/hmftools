package com.hartwig.hmftools.sage.context;

import java.util.Arrays;
import java.util.EnumSet;

import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.sage.read.IndexedBases;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class RefSequence {

    private static final EnumSet<CigarOperator> EXTEND_START = EnumSet.of(CigarOperator.I, CigarOperator.S);

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
        int additionalBases = 0;
        int cigarLength = 0;

        for (final CigarElement cigarElement : record.getCigar().getCigarElements()) {
            if (EXTEND_START.contains(cigarElement.getOperator())) {
                additionalBases += cigarElement.getLength();
            }
            cigarLength += cigarElement.getLength();
        }

        int refPositionStart = Math.max(1, alignmentStart - additionalBases);

        int alignmentEnd = Math.max(alignmentStart + cigarLength, alignmentStart + record.getReadLength());
        int alignmentStartIndex = alignmentStart - refPositionStart;

        return new IndexedBases(alignmentStart, alignmentStartIndex, alignment(refPositionStart, Math.min(alignmentEnd, sequenceLength)));
    }

}
