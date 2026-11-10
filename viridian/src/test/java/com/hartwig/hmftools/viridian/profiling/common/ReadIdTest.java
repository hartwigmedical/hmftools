package com.hartwig.hmftools.viridian.profiling.common;

import static com.hartwig.hmftools.viridian.ViridianTestUtils.samFlags;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;

public class ReadIdTest
{
    @Test
    public void testFromTakesMateFromPairFlags()
    {
        assertEquals(new ReadId("r1", 1), ReadId.from(record(samFlags(SAMFlag.READ_PAIRED, SAMFlag.FIRST_OF_PAIR))));
        assertEquals(new ReadId("r1", 2), ReadId.from(record(samFlags(SAMFlag.READ_PAIRED, SAMFlag.SECOND_OF_PAIR))));
        assertEquals(new ReadId("r1", 0), ReadId.from(record(samFlags())));
    }

    @Test
    public void testParseRoundTripsToString()
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
