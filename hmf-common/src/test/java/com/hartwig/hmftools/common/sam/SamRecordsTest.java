package com.hartwig.hmftools.common.sam;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SamRecordsTest {

    @Test
    public void testInsertedBasesAfterPosition() {
        assertEquals(0, SamRecords.basesInsertedAfterPosition(100, buildSamRecord(98, "6M", "GATACA")));
        assertEquals(1, SamRecords.basesInsertedAfterPosition(100, buildSamRecord(98, "3M1I3M", "GATTACA")));
        assertEquals(2, SamRecords.basesInsertedAfterPosition(100, buildSamRecord(98, "3M2I3M", "GATTTACA")));
        assertEquals(3, SamRecords.basesInsertedAfterPosition(100, buildSamRecord(98, "3M3I", "GATTTT")));

        assertEquals(0, SamRecords.basesInsertedAfterPosition(100, buildSamRecord(98, "6M", "GATACA")));
        assertEquals(0, SamRecords.basesInsertedAfterPosition(100, buildSamRecord(98, "3M1D2M", "GATCA")));
        assertEquals(0, SamRecords.basesInsertedAfterPosition(100, buildSamRecord(98, "3M2D1M", "GATA")));
        assertEquals(0, SamRecords.basesInsertedAfterPosition(100, buildSamRecord(98, "3M3D", "GAT")));
    }

    @Test
    public void testDeletedBasesAfterPosition() {
        assertEquals(0, SamRecords.basesDeletedAfterPosition(100, buildSamRecord(98, "6M", "GATACA")));
        assertEquals(0, SamRecords.basesDeletedAfterPosition(100, buildSamRecord(98, "3M1I3M", "GATTACA")));
        assertEquals(0, SamRecords.basesDeletedAfterPosition(100, buildSamRecord(98, "3M2I3M", "GATTTACA")));
        assertEquals(0, SamRecords.basesDeletedAfterPosition(100, buildSamRecord(98, "3M3I", "GATTTT")));

        assertEquals(0, SamRecords.basesDeletedAfterPosition(100, buildSamRecord(98, "6M", "GATACA")));
        assertEquals(1, SamRecords.basesDeletedAfterPosition(100, buildSamRecord(98, "3M1D2M", "GATCA")));
        assertEquals(2, SamRecords.basesDeletedAfterPosition(100, buildSamRecord(98, "3M2D1M", "GATA")));
        assertEquals(3, SamRecords.basesDeletedAfterPosition(100, buildSamRecord(98, "3M3D", "GAT")));
    }

    @NotNull
    public static SAMRecord buildSamRecord(final int alignmentStart, @NotNull final String cigar, @NotNull final String readString) {
        final StringBuilder qualityString = new StringBuilder();
        for (int i = 0; i < readString.length(); i++) {
            qualityString.append("A");
        }

        return SamRecordsTest.buildSamRecord(alignmentStart, cigar, readString, qualityString.toString());
    }

    @NotNull
    public static SAMRecord buildSamRecord(final int alignmentStart, @NotNull final String cigar, @NotNull final String readString,
            @NotNull final String qualities) {
        final SAMRecord record = new SAMRecord(null);
        record.setAlignmentStart(alignmentStart);
        record.setCigarString(cigar);
        record.setReadString(readString);
        record.setReadNegativeStrandFlag(false);
        record.setBaseQualityString(qualities);
        record.setMappingQuality(20);
        record.setDuplicateReadFlag(false);
        record.setReadUnmappedFlag(false);
        return record;
    }

}
