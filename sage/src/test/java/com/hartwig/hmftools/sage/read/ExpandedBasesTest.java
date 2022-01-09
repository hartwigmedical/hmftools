package com.hartwig.hmftools.sage.read;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.sage.evidence.ReadContextCounterTest;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ExpandedBasesTest
{

    private final ExpandedBasesFactory factory = new ExpandedBasesFactory(2, 3);

    @Test
    public void testNoSkipping()
    {
        String cigar = "3S2M1I1D2M1S";
        String read = "SSSMMIMMS";
        final SAMRecord samRecord = ReadContextCounterTest.buildSamRecord(100, cigar, read, read);

        for(int i = 0; i < read.length(); i++)
        {
            final IndexedBases victim = factory.expand(100, i, samRecord);
            assertExpand(i, read, victim);
        }
    }

    @Test
    public void testShortSkip()
    {
        String cigar = "2M1N2M";
        String read = "MMMM";
        final SAMRecord samRecord = ReadContextCounterTest.buildSamRecord(100, cigar, read, read);

        for(int i = 0; i < read.length(); i++)
        {
            final IndexedBases victim = factory.expand(100, i, samRecord);
            assertExpand(i, read, victim);
        }
    }

    @Test
    public void testOneSkip()
    {
        String cigar = "2M100N2M";
        String expectedRead = "MM...MM";
        String read = expectedRead.replaceAll("\\.", "");
        final SAMRecord samRecord = ReadContextCounterTest.buildSamRecord(100, cigar, read, read);

        assertExpand(0, expectedRead, factory.expand(100, 0, samRecord));
        assertExpand(1, expectedRead, factory.expand(100, 1, samRecord));
        assertExpand(5, expectedRead, factory.expand(100, 2, samRecord));
        assertExpand(6, expectedRead, factory.expand(100, 3, samRecord));
    }

    @Test
    public void testTwoSkip()
    {
        String cigar = "2M2456N2M100N3M";
        String expectedRead = "MM...MM...MMM";
        String read = expectedRead.replaceAll("\\.", "");
        final SAMRecord samRecord = ReadContextCounterTest.buildSamRecord(100, cigar, read, read);

        assertExpand(0, expectedRead, factory.expand(100, 0, samRecord));
        assertExpand(1, expectedRead, factory.expand(100, 1, samRecord));
        assertExpand(5, expectedRead, factory.expand(100, 2, samRecord));
        assertExpand(6, expectedRead, factory.expand(100, 3, samRecord));
        assertExpand(10, expectedRead, factory.expand(100, 4, samRecord));
        assertExpand(11, expectedRead, factory.expand(100, 5, samRecord));
        assertExpand(12, expectedRead, factory.expand(100, 6, samRecord));
    }

    @Test
    public void testThreeSkip()
    {
        String cigar = "2M100N1M100N3M100N4M";
        String expectedRead = "MM...M...MMM...MMMM";
        String read = expectedRead.replaceAll("\\.", "");
        final SAMRecord samRecord = ReadContextCounterTest.buildSamRecord(100, cigar, read, read);

        assertExpand(0, expectedRead, factory.expand(100, 0, samRecord));
        assertExpand(1, expectedRead, factory.expand(100, 1, samRecord));
        assertExpand(5, expectedRead, factory.expand(100, 2, samRecord));
        assertExpand(9, expectedRead, factory.expand(100, 3, samRecord));
        assertExpand(10, expectedRead, factory.expand(100, 4, samRecord));
        assertExpand(11, expectedRead, factory.expand(100, 5, samRecord));
        assertExpand(15, expectedRead, factory.expand(100, 6, samRecord));
        assertExpand(16, expectedRead, factory.expand(100, 7, samRecord));
        assertExpand(17, expectedRead, factory.expand(100, 8, samRecord));
    }

    private static void assertExpand(int expectedIndex, String expectedRead, IndexedBases victim)
    {
        assertEquals(expectedIndex, victim.Index);
        assertEquals(expectedRead, new String(victim.Bases));
    }

}
