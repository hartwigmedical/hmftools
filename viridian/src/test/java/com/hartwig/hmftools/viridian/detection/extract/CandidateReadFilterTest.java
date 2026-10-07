package com.hartwig.hmftools.viridian.detection.extract;

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

    private static final CandidateReadFilter FILTER =
            new CandidateReadFilter(CandidateReadSource.ALL, 20, Set.of("chrEBV"), null);

    private static final CandidateReadFilter FULLY_UNMAPPED_AND_DECOY_FILTER =
            new CandidateReadFilter(CandidateReadSource.FULLY_UNMAPPED_AND_VIRUS, 20, Set.of("chrEBV"), null);

    @Test
    public void testIsCandidateContigAllSourceAcceptsEveryContig()
    {
        assertTrue(FILTER.isCandidateContig("chr1"));
        assertTrue(FILTER.isCandidateContig("chrEBV"));
    }

    @Test
    public void testIsCandidateContigFullyUnmappedAndDecoySourceAcceptsOnlyDecoy()
    {
        assertFalse(FULLY_UNMAPPED_AND_DECOY_FILTER.isCandidateContig("chr1"));
        assertTrue(FULLY_UNMAPPED_AND_DECOY_FILTER.isCandidateContig("chrEBV"));
    }

    @Test
    public void testIsCandidateRecordFullyUnmappedAndDecoySourceRejectsHostContigCandidates()
    {
        List<SAMRecord> hostContigCandidates = List.of(
                read(0, "chr1", "30S70M"),
                read(73, "chr1", "100M"),
                readWithMate(0x1 | 0x40, "chr1", "100M", "chrEBV", 1000),
                readWithMate(0x1 | 0x4 | 0x80, "chr1", "*", "chr1", 100));

        for(SAMRecord record : hostContigCandidates)
        {
            assertTrue(FILTER.isCandidateRecord(record));
            assertFalse(FULLY_UNMAPPED_AND_DECOY_FILTER.isCandidateRecord(record));
        }
    }

    @Test
    public void testIsCandidateRecordFullyUnmappedAndDecoySourceAcceptsFullyUnmappedAndDecoyReads()
    {
        assertTrue(FULLY_UNMAPPED_AND_DECOY_FILTER.isCandidateRecord(read(4, "*", "*")));
        assertTrue(FULLY_UNMAPPED_AND_DECOY_FILTER.isCandidateRecord(read(0, "chrEBV", "100M")));
    }

    @Test
    public void testIsCandidateRecordRejectsDuplicateRead()
    {
        assertFalse(FILTER.isCandidateRecord(read(0x400, "chrEBV", "100M")));
    }

    @Test
    public void testIsCandidateRecordRejectsSecondaryAlignment()
    {
        assertFalse(FILTER.isCandidateRecord(read(0x100, "chrEBV", "100M")));
    }

    @Test
    public void testIsCandidateRecordRejectsSupplementaryAlignment()
    {
        assertFalse(FILTER.isCandidateRecord(read(0x800, "chrEBV", "100M")));
    }

    @Test
    public void testIsCandidateRecordAcceptsFullyUnmappedRead()
    {
        assertTrue(FILTER.isCandidateRecord(read(4, "*", "*")));
    }

    @Test
    public void testIsCandidateRecordRejectsReduxUnmappedFullyUnmappedRead()
    {
        SAMRecord record = read(4, "*", "*");
        record.setAttribute("UM", "chr1:100");
        assertFalse(FILTER.isCandidateRecord(record));
    }

    @Test
    public void testIsCandidateRecordAcceptsReduxUnmappedReadWithMateOnDecoy()
    {
        // Placed on the virus decoy by its mapped mate, so the fragment touches the virus and the unmapped-elsewhere
        // tag must not drop it.
        SAMRecord record = read(4, "chrEBV", "*");
        record.setAttribute("UM", "chr2:100");
        assertTrue(FILTER.isCandidateRecord(record));
    }

    @Test
    public void testIsCandidateRecordRejectsReduxUnmappedReadWithMateOnHost()
    {
        SAMRecord record = readWithMate(0x1 | 0x4 | 0x40, "chr1", "*", "chr1", 100);
        assertTrue(FILTER.isCandidateRecord(record));

        record.setAttribute("UM", "chr1:100");
        assertFalse(FILTER.isCandidateRecord(record));
    }

    @Test
    public void testIsCandidateRecordAcceptsMappedReadWithUnmappedMate()
    {
        // paired (0x1) + mate unmapped (0x8) + first of pair (0x40)
        assertTrue(FILTER.isCandidateRecord(read(73, "chr1", "100M")));
    }

    @Test
    public void testIsCandidateRecordAcceptsLongSoftClipEitherSide()
    {
        assertTrue(FILTER.isCandidateRecord(read(0, "chr1", "20S80M")));
        assertTrue(FILTER.isCandidateRecord(read(0, "chr1", "80M20S")));
    }

    @Test
    public void testIsCandidateRecordRejectsSoftClipBelowThreshold()
    {
        assertFalse(FILTER.isCandidateRecord(read(0, "chr1", "19S81M")));
    }

    @Test
    public void testIsCandidateRecordAcceptsDecoyMappedRead()
    {
        assertTrue(FILTER.isCandidateRecord(read(0, "chrEBV", "100M")));
    }

    @Test
    public void testIsCandidateRecordAcceptsHostReadWithMateOnDecoy()
    {
        // No significant clip, but the mate on the virus decoy anchors a host<->virus fragment.
        assertTrue(FILTER.isCandidateRecord(readWithMate(0x1 | 0x40, "chr1", "100M", "chrEBV", 1000)));
    }

    @Test
    public void testIsCandidateRecordRejectsHostReadWithMateOnHost()
    {
        assertFalse(FILTER.isCandidateRecord(readWithMate(0x1 | 0x40, "chr1", "100M", "chr1", 1000)));
    }

    @Test
    public void testIsCandidateRecordRejectsPlainMappedRead()
    {
        assertFalse(FILTER.isCandidateRecord(read(0, "chr1", "100M")));
    }

    @Test
    public void testIsCandidateRecordAcceptsClippedReadWithoutSupplementary()
    {
        // The clipped bases were not placed in the host, so the clip may mark a viral junction.
        assertTrue(FILTER.isCandidateRecord(read(0, "chr1", "30S70M")));
    }

    @Test
    public void testIsCandidateRecordRejectsClippedReadWithHostSupplementary()
    {
        // Clipped bases align elsewhere in the host, so the clip is not viral evidence.
        assertFalse(FILTER.isCandidateRecord(read(0, "chr1", "30S70M", "chr1,200,+,70M30S,60,0;")));
    }

    @Test
    public void testIsCandidateRecordAcceptsClippedReadWithViralSupplementary()
    {
        assertTrue(FILTER.isCandidateRecord(read(0, "chr1", "30S70M", "chrEBV,200,+,70M30S,60,0;")));
    }

    @Test
    public void testIsCandidateRecordRejectsClippedReadWithViralAndHostSupplementary()
    {
        assertFalse(FILTER.isCandidateRecord(read(0, "chr1", "30S70M", "chrEBV,200,+,70M30S,60,0;chr1,900,+,70M30S,60,0;")));
    }

    @Test
    public void testIsCandidateRecordAcceptsHostSupplementaryClipWithUnmappedMate()
    {
        // The supplementary check only gates the clip rule; the unmapped-mate rule still applies (0x1|0x8|0x40 = 73).
        assertTrue(FILTER.isCandidateRecord(read(73, "chr1", "30S70M", "chr1,200,+,70M30S,60,0;")));
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
