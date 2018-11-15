package com.hartwig.hmftools.common.sam;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class SAMRecords {

    public static final int PHRED_OFFSET = 33;

    public static int getBaseQuality(final char quality) {
        return quality - PHRED_OFFSET;
    }

    public static int getBaseQuality(@NotNull final SAMRecord record, int position) {
        return getAvgBaseQuality(record, position, 1);
    }

    public static int getAvgBaseQuality(@NotNull final SAMRecord record, int position, int length) {
        assert (position > 0);

        int score = 0;
        final String baseQualities = record.getBaseQualityString();
        for (int index = position - 1; index < Math.min(position - 1 + length, baseQualities.length()); index++) {
            score += getBaseQuality(baseQualities.charAt(index));
        }
        return score / length;
    }

    public static boolean containsInsert(@NotNull final SAMRecord record, int position, @NotNull final String alt) {
        int recordIdxOfVariantStart = record.getReadPositionAtReferencePosition(position);

        int insertedBases = basesInsertedAfterPosition(record, position);
        return insertedBases == alt.length() - 1 && record.getReadString()
                .substring(recordIdxOfVariantStart - 1, recordIdxOfVariantStart - 1 + alt.length())
                .equals(alt);
    }

    public static boolean containsDelete(@NotNull final SAMRecord record, int position, @NotNull final String ref) {
        int deletedBases = basesDeletedAfterPosition(record, position);
        return deletedBases == ref.length() - 1;
    }

    public static int basesInsertedAfterPosition(@NotNull final SAMRecord record, int position) {
        int startReadPosition = record.getReadPositionAtReferencePosition(position);
        assert (startReadPosition != 0);
        int nextReadPosition = record.getReadPositionAtReferencePosition(position + 1);

        return nextReadPosition == 0 && record.getAlignmentEnd() == position
                ? record.getReadString().length() - startReadPosition
                : Math.max(0, nextReadPosition - startReadPosition - 1);
    }

    public static int basesDeletedAfterPosition(@NotNull final SAMRecord record, int position) {
        int startReadPosition = record.getReadPositionAtReferencePosition(position);
        assert (startReadPosition != 0);

        int nextReferencePosition = record.getReferencePositionAtReadPosition(startReadPosition + 1);
        return nextReferencePosition == 0 && startReadPosition == record.getReadLength()
                ? record.getAlignmentEnd() - position
                : Math.max(0, nextReferencePosition - position - 1);
    }
}
