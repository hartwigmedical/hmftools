package com.hartwig.hmftools.common.utils.sam;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SAMRecordsTest {

    @Test
    public void testInsertedBasesAfterPosition() {
        assertEquals(0, SAMRecords.basesInsertedAfterPosition(buildSamRecord(98, "6M", "GATACA"), 100));
        assertEquals(1, SAMRecords.basesInsertedAfterPosition(buildSamRecord(98, "3M1I3M", "GATTACA"), 100));
        assertEquals(2, SAMRecords.basesInsertedAfterPosition(buildSamRecord(98, "3M2I3M", "GATTTACA"), 100));
        assertEquals(3, SAMRecords.basesInsertedAfterPosition(buildSamRecord(98, "3M3I", "GATTTT"), 100));

        assertEquals(0, SAMRecords.basesInsertedAfterPosition(buildSamRecord(98, "6M", "GATACA"), 100));
        assertEquals(0, SAMRecords.basesInsertedAfterPosition(buildSamRecord(98, "3M1D2M", "GATCA"), 100));
        assertEquals(0, SAMRecords.basesInsertedAfterPosition(buildSamRecord(98, "3M2D1M", "GATA"), 100));
        assertEquals(0, SAMRecords.basesInsertedAfterPosition(buildSamRecord(98, "3M3D", "GAT"), 100));
    }

    @Test
    public void testDeletedBasesAfterPosition() {
        assertEquals(0, SAMRecords.basesDeletedAfterPosition(buildSamRecord(98, "6M", "GATACA"), 100));
        assertEquals(0, SAMRecords.basesDeletedAfterPosition(buildSamRecord(98, "3M1I3M", "GATTACA"), 100));
        assertEquals(0, SAMRecords.basesDeletedAfterPosition(buildSamRecord(98, "3M2I3M", "GATTTACA"), 100));
        assertEquals(0, SAMRecords.basesDeletedAfterPosition(buildSamRecord(98, "3M3I", "GATTTT"), 100));

        assertEquals(0, SAMRecords.basesDeletedAfterPosition(buildSamRecord(98, "6M", "GATACA"), 100));
        assertEquals(1, SAMRecords.basesDeletedAfterPosition(buildSamRecord(98, "3M1D2M", "GATCA"), 100));
        assertEquals(2, SAMRecords.basesDeletedAfterPosition(buildSamRecord(98, "3M2D1M", "GATA"), 100));
        assertEquals(3, SAMRecords.basesDeletedAfterPosition(buildSamRecord(98, "3M3D", "GAT"), 100));
    }

    @NotNull
    public static SAMRecord buildSamRecord(final int alignmentStart, @NotNull final String cigar, @NotNull final String readString) {
        final StringBuilder qualityString = new StringBuilder();
        for (int i = 0; i < readString.length(); i++) {
            qualityString.append("A");
        }

        return SAMRecordsTest.buildSamRecord(alignmentStart, cigar, readString, qualityString.toString());
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
