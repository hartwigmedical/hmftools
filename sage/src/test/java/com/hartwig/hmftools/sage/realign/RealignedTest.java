package com.hartwig.hmftools.sage.realign;

import static com.hartwig.hmftools.sage.realign.Realigned.realigned;
import static com.hartwig.hmftools.sage.realign.RealignedType.EXACT;
import static com.hartwig.hmftools.sage.realign.RealignedType.LENGTHENED;
import static com.hartwig.hmftools.sage.realign.RealignedType.NONE;
import static com.hartwig.hmftools.sage.realign.RealignedType.SHORTENED;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class RealignedTest
{

    @Test
    public void testRealignedTooShort()
    {
        String sequence = "GAGAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTTT";
        String truncatedAtEnd = "GAGAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTT";
        String truncatedAtStart = "GAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTTT";
        int startIndex = 0;
        int endIndex = sequence.length() - 1;

        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex, sequence.getBytes(), 0));
        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), startIndex, truncatedAtEnd.getBytes(), 0));
        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), -2, truncatedAtStart.getBytes(), 0));
    }

    @Test
    public void testRealigned()
    {
        String sequence = "GAGAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTTTATTATTG";
        int startIndex = 3;
        int endIndex = startIndex + 55;

        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex, sequence.getBytes(), 0));

        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), startIndex - 1, sequence.getBytes(), 0));
        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), startIndex + 1, sequence.getBytes(), 0));

        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), startIndex - 2, sequence.getBytes(), 1));
        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex - 1, sequence.getBytes(), 1));
        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex + 1, sequence.getBytes(), 1));
        assertRealigned(NONE, realigned(startIndex, endIndex, sequence.getBytes(), startIndex + 2, sequence.getBytes(), 1));

        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex - 2, sequence.getBytes(), 2));
        assertRealigned(EXACT, realigned(startIndex, endIndex, sequence.getBytes(), startIndex + 2, sequence.getBytes(), 2));
    }

    @Test
    public void testPolyA()
    {
        String shorter = "GATCAAAAAAAAAGATC";
        String ref = "GATCAAAAAAAAAAGATC";
        String longer = "GATCAAAAAAAAAAAGATC";

        int startIndex = 0;
        int endIndex = ref.length() - 1;

        assertRealigned(EXACT, realigned(startIndex, endIndex, ref.getBytes(), startIndex, ref.getBytes(), 10));
        assertRealigned(SHORTENED, 9, realigned(startIndex, endIndex, ref.getBytes(), startIndex, shorter.getBytes(), 10));
        assertRealigned(LENGTHENED, 11, realigned(startIndex, endIndex, ref.getBytes(), startIndex, longer.getBytes(), 10));
    }

    @Test
    public void testDiNucleotideRepeat()
    {
        String shorter = "GATCATATATATGATC";
        String ref = "GATCATATATATATGATC";
        String longer = "GATCATATATATATATGATC";

        int startIndex = 0;
        int endIndex = ref.length() - 1;

        assertRealigned(EXACT, realigned(startIndex, endIndex, ref.getBytes(), startIndex, ref.getBytes(), 10));
        assertRealigned(SHORTENED, 4, realigned(startIndex, endIndex, ref.getBytes(), startIndex, shorter.getBytes(), 10));
        assertRealigned(LENGTHENED, 6, realigned(startIndex, endIndex, ref.getBytes(), startIndex, longer.getBytes(), 10));
    }

    private static void assertRealigned(RealignedType expectedType, int expectedCount, RealignedContext context)
    {
        assertEquals(expectedCount, context.RepeatCount);
        assertEquals(expectedType, context.Type);
    }

    private static void assertRealigned(RealignedType expectedType, RealignedContext context)
    {
        assertEquals(expectedType, context.Type);
    }
}
