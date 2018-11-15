package com.hartwig.hmftools.common.sam;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;

public class SamRecords {

    public static int basesInsertedAfterPosition(int position, @NotNull final SAMRecord record) {
        int startReadPosition = record.getReadPositionAtReferencePosition(position);
        assert (startReadPosition != 0);
        int nextReadPosition = record.getReadPositionAtReferencePosition(position + 1);

        return nextReadPosition == 0 && record.getAlignmentEnd() == position
                ? record.getReadString().length() - startReadPosition
                : Math.max(0, nextReadPosition - startReadPosition - 1);
    }

    public static int basesDeletedAfterPosition(int position, @NotNull final SAMRecord record) {
        int startReadPosition = record.getReadPositionAtReferencePosition(position);
        assert (startReadPosition != 0);

        int nextReferencePosition = record.getReferencePositionAtReadPosition(startReadPosition + 1);
        return nextReferencePosition == 0 && startReadPosition == record.getReadLength()
                ? record.getAlignmentEnd() - position
                : Math.max(0, nextReferencePosition - position - 1);
    }
}
