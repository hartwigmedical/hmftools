package com.hartwig.hmftools.sage.common;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.sage.common.NumberEvents;
import com.hartwig.hmftools.sage.common.RefSequence;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequence;

public class NumberEventsTest
{
    @Test
    public void testRawNumberEvents()
    {
        int refStart = 4000;
        SAMRecord record = buildSamRecord(refStart + 4, "7M", "GATTACA");
        assertEquals(0, NumberEvents.rawNM(record, new RefSequence(new ReferenceSequence("1",
                refStart - 1, "AAAAGATTACAGGGGG".getBytes()))));
        assertEquals(1, NumberEvents.rawNM(record, new RefSequence(new ReferenceSequence("1",
                refStart - 1, "AAAAGATTAAAGGGG".getBytes()))));
        assertEquals(1, NumberEvents.rawNM(record, new RefSequence(new ReferenceSequence("1",
                refStart - 1, "AAAAGATTACCGGGG".getBytes()))));
        assertEquals(2, NumberEvents.rawNM(record, new RefSequence(new ReferenceSequence("1",
                refStart - 1, "AAAAGATTAGGGGGG".getBytes()))));
    }

    private static SAMRecord buildSamRecord(final int alignmentStart, @NotNull final String cigar, @NotNull final String readString)
    {
        final StringBuilder qualityString = new StringBuilder();
        for(int i = 0; i < readString.length(); i++)
        {
            qualityString.append("A");
        }

        return buildSamRecord(alignmentStart, cigar, readString, qualityString.toString());
    }

    @NotNull
    public static SAMRecord buildSamRecord(final int alignmentStart, @NotNull final String cigar, @NotNull final String readString,
            @NotNull final String qualities)
    {
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
