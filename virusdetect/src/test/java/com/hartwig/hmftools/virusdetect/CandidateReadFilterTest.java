package com.hartwig.hmftools.virusdetect;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import org.junit.Test;

public class CandidateReadFilterTest
{
    private static final SAMSequenceDictionary DICT = new SAMSequenceDictionary(List.of(
            new SAMSequenceRecord("chr1", 10000),
            new SAMSequenceRecord("chrEBV", 171823)));

    private static final String BASES = "A".repeat(100);

    private static final CandidateReadFilter FILTER = new CandidateReadFilter(20, Set.of("chrEBV"));

    @Test
    public void testUnmappedReadIsCandidate()
    {
        assertTrue(FILTER.isCandidate(read(4, "*", "*")));
    }

    @Test
    public void testMappedReadWithUnmappedMateIsCandidate()
    {
        // paired (0x1) + mate unmapped (0x8) + first of pair (0x40)
        assertTrue(FILTER.isCandidate(read(73, "chr1", "100M")));
    }

    @Test
    public void testLongSoftClipIsCandidateEitherSide()
    {
        assertTrue(FILTER.isCandidate(read(0, "chr1", "20S80M")));
        assertTrue(FILTER.isCandidate(read(0, "chr1", "80M20S")));
    }

    @Test
    public void testSoftClipBelowThresholdIsNotCandidate()
    {
        assertFalse(FILTER.isCandidate(read(0, "chr1", "19S81M")));
    }

    @Test
    public void testDecoyMappedReadIsCandidate()
    {
        assertTrue(FILTER.isCandidate(read(0, "chrEBV", "100M")));
    }

    @Test
    public void testPlainMappedReadIsNotCandidate()
    {
        assertFalse(FILTER.isCandidate(read(0, "chr1", "100M")));
    }

    private static SAMRecord read(int flags, String contig, String cigar)
    {
        String position = contig.equals("*") ? "0" : "100";
        String line = String.join("\t", "read", String.valueOf(flags), contig, position, "0", cigar, "*", "0", "0", BASES, "*");
        return SamRecordTestUtils.parseSamString(line, DICT);
    }
}
