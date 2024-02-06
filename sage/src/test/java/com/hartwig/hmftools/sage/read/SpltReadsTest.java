package com.hartwig.hmftools.sage.read;

import static com.hartwig.hmftools.sage.common.TestUtils.buildSamRecord;
import static com.hartwig.hmftools.sage.read.SplitReadUtils.MAX_SKIPPED_REFERENCE_REGIONS;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.sage.evidence.ReadIndexBases;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SpltReadsTest
{
    private final String SKIP_WILDCARDS;

    public SpltReadsTest()
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
        final SAMRecord samRecord = buildSamRecord(100, cigar, read, read);

        for(int i = 0; i < read.length(); i++)
        {
            final ReadIndexBases readIndexBases = SplitReadUtils.expandSplitRead(i, samRecord);
            assertExpand(i, read, readIndexBases);
        }
    }

    @Test
    public void testShortSkip()
    {
        String cigar = "2M1N2M";
        String read = "MMMM";
        final SAMRecord samRecord = buildSamRecord(100, cigar, read, read);

        for(int i = 0; i < read.length(); i++)
        {
            final ReadIndexBases readIndexBases = SplitReadUtils.expandSplitRead(i, samRecord);
            assertExpand(i, read, readIndexBases);
        }
    }

    @Test
    public void testOneSkip()
    {
        String cigar = "2M100N2M";
        String expectedRead = "MM" + SKIP_WILDCARDS + "MM";
        String read = expectedRead.replaceAll("\\.", "");
        final SAMRecord samRecord = buildSamRecord(100, cigar, read, read);

        assertExpand(0, expectedRead, SplitReadUtils.expandSplitRead(0, samRecord));
        assertExpand(1, expectedRead, SplitReadUtils.expandSplitRead(1, samRecord));
        assertExpand(2 + MAX_SKIPPED_REFERENCE_REGIONS, expectedRead, SplitReadUtils.expandSplitRead(2, samRecord));
        assertExpand(3 + MAX_SKIPPED_REFERENCE_REGIONS, expectedRead, SplitReadUtils.expandSplitRead(3, samRecord));
    }

    @Test
    public void testTwoSkip()
    {
        String cigar = "2M2456N2M100N3M";
        String expectedRead = "MM" + SKIP_WILDCARDS + "MM" + SKIP_WILDCARDS + "MMM";
        String read = expectedRead.replaceAll("\\.", "");
        final SAMRecord samRecord = buildSamRecord(100, cigar, read, read);

        assertExpand(0, expectedRead, SplitReadUtils.expandSplitRead(0, samRecord));
        assertExpand(1, expectedRead, SplitReadUtils.expandSplitRead(1, samRecord));
        assertExpand(2 + MAX_SKIPPED_REFERENCE_REGIONS, expectedRead, SplitReadUtils.expandSplitRead(2, samRecord));
        assertExpand(3 + MAX_SKIPPED_REFERENCE_REGIONS, expectedRead, SplitReadUtils.expandSplitRead(3, samRecord));
        assertExpand(4 + MAX_SKIPPED_REFERENCE_REGIONS * 2, expectedRead, SplitReadUtils.expandSplitRead(4, samRecord));
        assertExpand(5 + MAX_SKIPPED_REFERENCE_REGIONS * 2, expectedRead, SplitReadUtils.expandSplitRead(5, samRecord));
    }

    private static void assertExpand(int expectedIndex, final String expectedRead, final ReadIndexBases readIndexBases)
    {
        assertEquals(expectedIndex, readIndexBases.Index);
        assertEquals(expectedRead, new String(readIndexBases.Bases));
    }

}
