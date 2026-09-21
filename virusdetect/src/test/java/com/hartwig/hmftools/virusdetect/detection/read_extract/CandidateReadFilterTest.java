package com.hartwig.hmftools.virusdetect.detection.read_extract;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.jetbrains.annotations.Nullable;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class CandidateReadFilterTest
{
    private static final SAMSequenceDictionary DICT = new SAMSequenceDictionary(List.of(
            new SAMSequenceRecord("chr1", 10000),
            new SAMSequenceRecord("chrEBV", 171823)));

    private static final String BASES = "A".repeat(100);

    private static final CandidateReadFilter FILTER = new CandidateReadFilter(20, Set.of("chrEBV"));

    @Test
    public void testDuplicateReadIsNotCandidate()
    {
        assertFalse(FILTER.isCandidate(read(0x400, "chrEBV", "100M")));
    }

    @Test
    public void testSecondaryAlignmentIsNotCandidate()
    {
        assertFalse(FILTER.isCandidate(read(0x100, "chrEBV", "100M")));
    }

    @Test
    public void testSupplementaryAlignmentIsNotCandidate()
    {
        assertFalse(FILTER.isCandidate(read(0x800, "chrEBV", "100M")));
    }

    @Test
    public void testUnmappedReadIsCandidate()
    {
        assertTrue(FILTER.isCandidate(read(4, "*", "*")));
    }

    @Test
    public void testReduxUnmappedReadIsNotCandidate()
    {
        SAMRecord record = read(4, "*", "*");
        record.setAttribute("UM", "chr1:100");
        assertFalse(FILTER.isCandidate(record));
    }

    @Test
    public void testReduxUnmappedReadPlacedOnDecoyIsCandidate()
    {
        // Placed on the viral decoy by its mapped mate, so the fragment touches the virus and the unmapped-elsewhere
        // tag must not drop it.
        SAMRecord record = read(4, "chrEBV", "*");
        record.setAttribute("UM", "chr2:100");
        assertTrue(FILTER.isCandidate(record));
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
    public void testHostReadWithMateOnDecoyIsCandidate()
    {
        // No significant clip, but the mate on the viral decoy anchors a host<->virus fragment.
        assertTrue(FILTER.isCandidate(readWithMate(0x1 | 0x40, "chr1", "100M", "chrEBV", 1000)));
    }

    @Test
    public void testHostReadWithMateOnHostIsNotCandidate()
    {
        assertFalse(FILTER.isCandidate(readWithMate(0x1 | 0x40, "chr1", "100M", "chr1", 1000)));
    }

    @Test
    public void testPlainMappedReadIsNotCandidate()
    {
        assertFalse(FILTER.isCandidate(read(0, "chr1", "100M")));
    }

    @Test
    public void testClippedReadWithoutSupplementaryIsCandidate()
    {
        // The clipped bases were not placed in the host, so the clip may mark a viral junction.
        assertTrue(FILTER.isCandidate(read(0, "chr1", "30S70M")));
    }

    @Test
    public void testClippedReadWithHostSupplementaryIsNotCandidate()
    {
        // Clipped bases align elsewhere in the host, so the clip is not viral evidence.
        assertFalse(FILTER.isCandidate(read(0, "chr1", "30S70M", "chr1,200,+,70M30S,60,0;")));
    }

    @Test
    public void testClippedReadWithViralSupplementaryIsCandidate()
    {
        assertTrue(FILTER.isCandidate(read(0, "chr1", "30S70M", "chrEBV,200,+,70M30S,60,0;")));
    }

    @Test
    public void testClippedReadWithViralAndHostSupplementaryIsNotCandidate()
    {
        assertFalse(FILTER.isCandidate(read(0, "chr1", "30S70M", "chrEBV,200,+,70M30S,60,0;chr1,900,+,70M30S,60,0;")));
    }

    @Test
    public void testHostSupplementaryClipWithUnmappedMateIsStillCandidate()
    {
        // The supplementary check only gates the clip rule; the unmapped-mate rule still applies (0x1|0x8|0x40 = 73).
        assertTrue(FILTER.isCandidate(read(73, "chr1", "30S70M", "chr1,200,+,70M30S,60,0;")));
    }

    private static SAMRecord read(int flags, String contig, String cigar)
    {
        return read(flags, contig, cigar, null);
    }

    private static SAMRecord readWithMate(int flags, String contig, String cigar, String mateContig, int matePosition)
    {
        String position = contig.equals("*") ? "0" : "100";
        String line = String.join(
                "\t", "read", String.valueOf(flags), contig, position, "0", cigar,
                mateContig, String.valueOf(matePosition), "0", BASES, "*");
        return SamRecordTestUtils.parseSamString(line, DICT);
    }

    private static SAMRecord read(int flags, String contig, String cigar, @Nullable String supplementaryTag)
    {
        String position = contig.equals("*") ? "0" : "100";
        String line = String.join("\t", "read", String.valueOf(flags), contig, position, "0", cigar, "*", "0", "0", BASES, "*");
        if(supplementaryTag != null)
        {
            line = line + "\tSA:Z:" + supplementaryTag;
        }
        return SamRecordTestUtils.parseSamString(line, DICT);
    }
}
