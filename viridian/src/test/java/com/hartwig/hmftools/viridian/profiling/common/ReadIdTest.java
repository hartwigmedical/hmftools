package com.hartwig.hmftools.viridian.profiling.common;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ReadIdTest
{
    @Test
    public void testFromRecordTakesMateFromPairFlags()
    {
        assertEquals(new ReadId("r1", 1), ReadId.from(record(0x1 | 0x40)));
        assertEquals(new ReadId("r1", 2), ReadId.from(record(0x1 | 0x80)));
        assertEquals(new ReadId("r1", 0), ReadId.from(record(0)));
    }

    @Test
    public void testStringRepresentationRoundTrips()
    {
        assertEquals("r1/1", new ReadId("r1", 1).toString());
        assertEquals("r1/2", new ReadId("r1", 2).toString());
        assertEquals("r1", new ReadId("r1", 0).toString());

        for(ReadId readId : new ReadId[] { new ReadId("r1", 1), new ReadId("r1", 2), new ReadId("r1", 0) })
        {
            assertEquals(readId, ReadId.parse(readId.toString()));
        }
    }

    private static SAMRecord record(int flags)
    {
        SAMRecord record = new SAMRecord(null);
        record.setReadName("r1");
        record.setFlags(flags);
        return record;
    }
}
