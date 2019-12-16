package com.hartwig.hmftools.sage.context;

import static org.junit.Assert.assertFalse;

import static junit.framework.TestCase.assertTrue;

import org.junit.Test;

public class RealignedTest {

    @Test
    public void testRealignedTooShort() {
        String sequence = "GAGAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTTT";
        String truncatedAtEnd = "GAGAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTT";
        String truncatedAtStart = "GAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTTT";
        int startIndex = 3;
        int endIndex = startIndex + 55;

        assertTrue(Realigned.realigned(startIndex, startIndex-2, endIndex, sequence.getBytes(), startIndex, sequence.getBytes(), 0));
        assertFalse(Realigned.realigned(startIndex, startIndex-2, endIndex, sequence.getBytes(), startIndex, truncatedAtEnd.getBytes(), 0));
        assertFalse(Realigned.realigned(startIndex, startIndex-2, endIndex, sequence.getBytes(), startIndex - 2, truncatedAtStart.getBytes(), 0));

    }

    @Test
    public void testRealigned() {
        String sequence = "GAGAGTGTGTGTGTGTGTCTGTGTGTATGTATATATATATATATATATATCACATTTTTATTATTG";
        int startIndex = 3;
        int endIndex = startIndex + 55;

        assertTrue(Realigned.realigned(startIndex, startIndex, endIndex, sequence.getBytes(), startIndex, sequence.getBytes(), 0));

        assertFalse(Realigned.realigned(startIndex, startIndex, endIndex, sequence.getBytes(), startIndex - 1, sequence.getBytes(), 0));
        assertFalse(Realigned.realigned(startIndex, startIndex, endIndex, sequence.getBytes(), startIndex + 1, sequence.getBytes(), 0));
        assertTrue(Realigned.realigned(startIndex, startIndex, endIndex, sequence.getBytes(), startIndex - 1, sequence.getBytes(), 1));
        assertTrue(Realigned.realigned(startIndex, startIndex, endIndex, sequence.getBytes(), startIndex + 1, sequence.getBytes(), 1));

        assertFalse(Realigned.realigned(startIndex, startIndex, endIndex, sequence.getBytes(), startIndex - 2, sequence.getBytes(), 1));
        assertFalse(Realigned.realigned(startIndex, startIndex, endIndex, sequence.getBytes(), startIndex + 2, sequence.getBytes(), 1));
        assertTrue(Realigned.realigned(startIndex, startIndex, endIndex, sequence.getBytes(), startIndex - 2, sequence.getBytes(), 2));
        assertTrue(Realigned.realigned(startIndex, startIndex, endIndex, sequence.getBytes(), startIndex + 2, sequence.getBytes(), 2));
    }

}
