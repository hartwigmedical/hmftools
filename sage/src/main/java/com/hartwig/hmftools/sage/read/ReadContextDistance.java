package com.hartwig.hmftools.sage.read;

import com.google.common.annotations.VisibleForTesting;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ReadContextDistance {

    private int distance;
    private final String difference;

    @VisibleForTesting
    ReadContextDistance(final boolean testing, final int buffer, final int readBasesAltIndex, @NotNull final SAMRecord record,
            byte[] refBases) {
        this(readBasesAltIndex - buffer, readBasesAltIndex + buffer, record, new IndexedBases(record.getAlignmentStart(), 0, refBases));
    }

    public ReadContextDistance(final int readBasesStartIndex, final int readBasesEndIndex, @NotNull final SAMRecord record, @NotNull final IndexedBases refBases) {

        final StringBuilder differenceBuilder = new StringBuilder();

        int refBasesCurrentIndex = refBases.index(record.getAlignmentStart());
        int readBasesCurrentIndex = 0;
        for (CigarElement cigarElement : record.getCigar()) {

            final CigarOperator operator = cigarElement.getOperator();

            int refBasesNextIndex = refBasesCurrentIndex;
            int readBasesNextIndex = readBasesCurrentIndex;

            // Past end
            if (readBasesCurrentIndex > readBasesEndIndex) {
                break;
            }

            if (operator.consumesReadBases()) {
                readBasesNextIndex += cigarElement.getLength();
            }

            if (operator.consumesReferenceBases()) {
                refBasesNextIndex += cigarElement.getLength();
            }

            boolean isOverlap = readBasesStartIndex <= readBasesNextIndex - 1;
            if (isOverlap) {

                int readBasesMinIndex = Math.max(readBasesStartIndex, readBasesCurrentIndex);
                int readBasesMaxIndex = Math.min(readBasesNextIndex - 1, readBasesEndIndex);
                int readBasesLength = readBasesMaxIndex - readBasesMinIndex + 1;

                int refBasesMinIndex = readBasesMinIndex - readBasesCurrentIndex + refBasesCurrentIndex;

                switch (cigarElement.getOperator()) {
                    case S:
                    case I:
                    case X:
                        distance++;
                        differenceBuilder.append(readBasesLength).append(cigarElement.getOperator().toString());
                        break;
                    case D:
                        distance++;
                        differenceBuilder.append(cigarElement.getLength()).append(cigarElement.getOperator().toString());
                        break;
                    case EQ:
                        differenceBuilder.append(readBasesLength).append(CigarOperator.EQ.toString());
                        break;
                    case M:
                        final String cigar = mCigar(readBasesLength, readBasesMinIndex, record.getReadBases(), refBasesMinIndex, refBases.bases());
                        differenceBuilder.append(cigar);
                        break;

                }
            }

            readBasesCurrentIndex = readBasesNextIndex;
            refBasesCurrentIndex = refBasesNextIndex;

        }

        difference = differenceBuilder.toString();
    }

    @NotNull
    private String mCigar(int length, int readBaseIndex, byte[] readBases, int refBaseIndex, byte[] refBases) {

        final StringBuilder differenceBuilder = new StringBuilder();

        char currentOperator = 'X';
        int currentOperatorCount = 0;

        for (int i = 0; i < length; i++) {

            byte readBase = readBases[readBaseIndex + i];
            byte refBase = refBases[refBaseIndex + i];

            char newOperator;
            if (readBase == refBase) {
                newOperator = 'M';
            } else {
                newOperator = 'X';
                distance++;
            }

            if (currentOperatorCount > 0 && newOperator != currentOperator) {
                differenceBuilder.append(currentOperatorCount).append(currentOperator);
                currentOperatorCount = 0;
            }

            currentOperator = newOperator;
            currentOperatorCount++;
        }

        if (currentOperatorCount > 0) {
            differenceBuilder.append(currentOperatorCount).append(currentOperator);
        }

        return differenceBuilder.toString();
    }

    @NotNull
    public String cigar() {
        return difference;
    }

    public int distance() {
        return distance;
    }
}
