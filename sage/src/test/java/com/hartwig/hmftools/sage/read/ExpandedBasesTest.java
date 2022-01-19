package com.hartwig.hmftools.sage.read;

import static com.hartwig.hmftools.sage.read.ExpandedBasesFactory.MAX_SKIPPED_REFERENCE_REGIONS;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.sage.common.IndexedBases;
import com.hartwig.hmftools.sage.evidence.ReadContextCounterTest;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class ExpandedBasesTest
{
    // private final ExpandedBasesFactory factory = new ExpandedBasesFactory(2, 3);

    private final String SKIP_WILDCARDS;

    public ExpandedBasesTest()
    {
        StringBuilder sb = new StringBuilder();

        for(int i = 0; i < MAX_SKIPPED_REFERENCE_REGIONS; ++i)
        {
            sb.append('.');
        }

        SKIP_WILDCARDS = sb.toString();
    }

    @Test
    public void testNoSkipping()
    {
        String cigar = "3S2M1I1D2M1S";
        String read = "SSSMMIMMS";
        final SAMRecord samRecord = ReadContextCounterTest.buildSamRecord(100, cigar, read, read);

        for(int i = 0; i < read.length(); i++)
        {
            final IndexedBases victim = ExpandedBasesFactory.expand(100, i, samRecord);
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
            final IndexedBases victim = ExpandedBasesFactory.expand(100, i, samRecord);
            assertExpand(i, read, victim);
        }
    }

    @Test
    public void testOneSkip()
    {
        String cigar = "2M100N2M";
        String expectedRead = "MM" + SKIP_WILDCARDS + "MM";
        String read = expectedRead.replaceAll("\\.", "");
        final SAMRecord samRecord = ReadContextCounterTest.buildSamRecord(100, cigar, read, read);

        assertExpand(0, expectedRead, ExpandedBasesFactory.expand(100, 0, samRecord));
        assertExpand(1, expectedRead, ExpandedBasesFactory.expand(100, 1, samRecord));
        assertExpand(2 + MAX_SKIPPED_REFERENCE_REGIONS, expectedRead, ExpandedBasesFactory.expand(100, 2, samRecord));
        assertExpand(3 + MAX_SKIPPED_REFERENCE_REGIONS, expectedRead, ExpandedBasesFactory.expand(100, 3, samRecord));
    }

    @Test
    public void testTwoSkip()
    {
        String cigar = "2M2456N2M100N3M";
        String expectedRead = "MM" + SKIP_WILDCARDS + "MM" + SKIP_WILDCARDS + "MMM";
        String read = expectedRead.replaceAll("\\.", "");
        final SAMRecord samRecord = ReadContextCounterTest.buildSamRecord(100, cigar, read, read);

        assertExpand(0, expectedRead, ExpandedBasesFactory.expand(100, 0, samRecord));
        assertExpand(1, expectedRead, ExpandedBasesFactory.expand(100, 1, samRecord));
        assertExpand(2 + MAX_SKIPPED_REFERENCE_REGIONS, expectedRead, ExpandedBasesFactory.expand(100, 2, samRecord));
        assertExpand(3 + MAX_SKIPPED_REFERENCE_REGIONS, expectedRead, ExpandedBasesFactory.expand(100, 3, samRecord));
        assertExpand(4 + MAX_SKIPPED_REFERENCE_REGIONS * 2, expectedRead, ExpandedBasesFactory.expand(100, 4, samRecord));
        assertExpand(5 + MAX_SKIPPED_REFERENCE_REGIONS * 2, expectedRead, ExpandedBasesFactory.expand(100, 5, samRecord));
    }

    private static void assertExpand(int expectedIndex, String expectedRead, IndexedBases victim)
    {
        assertEquals(expectedIndex, victim.Index);
        assertEquals(expectedRead, new String(victim.Bases));
    }

}
