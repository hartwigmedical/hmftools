package com.hartwig.hmftools.virusdetect;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.TextCigarCodec;

public class ContigStatsCalculatorTest
{
    @Rule
    public TemporaryFolder mTempDir = new TemporaryFolder();

    private static final double EPSILON = 1e-9;

    // v1 (length 20): read r1 covers 1-10 (score 10), read r2 covers 6-10 (score 8); a lower-scoring
    // secondary of r1 on v1 is deduped away. v2 (length 10): read r3 covers it fully.
    @Test
    public void testComputesPerContigDepthAndCoverage() throws IOException
    {
        ViralReference reference = reference();
        SAMFileHeader header = header();

        List<SAMRecord> records = List.of(
                aligned(header, "r1", 0, "v1", 1, "10M", 10),
                aligned(header, "r1", 0x100, "v1", 1, "10M", 5),   // lower-scoring duplicate on same contig
                aligned(header, "r2", 0, "v1", 6, "5M", 8),
                aligned(header, "r3", 0, "v2", 1, "10M", 9));

        Map<String, ContigStats> stats = new ContigStatsCalculator().compute(writeBam(header, records), reference);

        assertEquals(2, stats.size());

        ContigStats v1 = stats.get("v1");
        assertEquals(20, v1.contigLength());
        assertEquals(2, v1.readCount());              // r1 counted once despite two alignments
        assertEquals(10, v1.coveredBases());          // positions 1-10
        assertEquals(0, v1.minDepth());               // positions 11-20 uncovered
        assertEquals(2, v1.maxDepth());               // positions 6-10 covered by both reads
        assertEquals(15.0 / 20, v1.meanDepth(), EPSILON);
        assertEquals(9.0, v1.meanAlignerScore(), EPSILON);   // (10 + 8) / 2
        assertEquals(0.5, v1.coverageFraction(), EPSILON);

        ContigStats v2 = stats.get("v2");
        assertEquals(1, v2.readCount());
        assertEquals(10, v2.coveredBases());
        assertEquals(1, v2.minDepth());
        assertEquals(1, v2.maxDepth());
        assertEquals(1.0, v2.meanDepth(), EPSILON);
        assertEquals(9.0, v2.meanAlignerScore(), EPSILON);
        assertEquals(1.0, v2.coverageFraction(), EPSILON);
    }

    // r1 and r2 each align to both contigs, more closely to v1 (fewer mismatches). Both reads' votes and their
    // contested margins land on v1; v2, never a read's best, holds no margins.
    @Test
    public void testVotesAndMarginsAttributeStrainSupport() throws IOException
    {
        ViralReference reference = reference();
        SAMFileHeader header = header();

        List<SAMRecord> records = List.of(
                aligned(header, "r1", 0, "v1", 1, "10M", 10, 1),
                aligned(header, "r1", 0x100, "v2", 1, "10M", 6, 4),
                aligned(header, "r2", 0, "v1", 1, "10M", 10, 0),
                aligned(header, "r2", 0x100, "v2", 1, "10M", 5, 5));

        // Injected correct-base probability 0.5, so each extra divergent base halves a contig's weight (0.5^diff).
        Map<String, ContigStats> stats = new ContigStatsCalculator(0.5).compute(writeBam(header, records), reference);

        ContigStats v1 = stats.get("v1");
        ContigStats v2 = stats.get("v2");

        // r1: v1 weight 0.5^0, v2 0.5^3; r2: v1 0.5^0, v2 0.5^5. Votes per read sum to 1, so both contigs sum to 2.
        assertEquals(1.0 / (1 + Math.pow(0.5, 3)) + 1.0 / (1 + Math.pow(0.5, 5)), v1.voteReads(), EPSILON);
        assertEquals(2.0 - v1.voteReads(), v2.voteReads(), EPSILON);

        // v1 wins both reads; margins are runner-up minus best edit distance: r1 4-1=3, r2 5-0=5.
        assertEquals(2, v1.contestedReads());
        assertEquals(4.0, v1.marginMean(), EPSILON);
        assertEquals(3.0, v1.marginMedian(), EPSILON);
        assertEquals(5.0, v1.marginP90(), EPSILON);

        // v2 is never a read's best, so it holds no margins.
        assertEquals(0, v2.contestedReads());
        assertTrue(Double.isNaN(v2.marginMean()));
    }

    // A read matches v1 full-length with 5 mismatches, but v2 only after clipping 40 bases (1 mismatch in the rest).
    // Counting the clip against v2, v1 explains the read far better and wins, despite v2's lower NM.
    @Test
    public void testClippingCountsAgainstContigInRivalry() throws IOException
    {
        List<ViralContig> contigs = List.of(
                new ViralContig("v1", 100, "Virus 1", "Group 1"),
                new ViralContig("v2", 100, "Virus 1", "Group 1"));
        ViralReference reference = new ViralReference(contigs, new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("v1", 100), new SAMSequenceRecord("v2", 100))));

        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        header.addSequence(new SAMSequenceRecord("v1", 100));
        header.addSequence(new SAMSequenceRecord("v2", 100));

        List<SAMRecord> records = List.of(
                aligned(header, "r", 0, "v1", 1, "100M", 90, 5),
                aligned(header, "r", 0x100, "v2", 1, "40S60M", 60, 1));

        Map<String, ContigStats> stats = new ContigStatsCalculator(0.5).compute(writeBam(header, records), reference);

        assertEquals(1, stats.get("v1").contestedReads());
        assertEquals(0, stats.get("v2").contestedReads());
        // margin = v2 divergence (1 mismatch + 40 clipped) - v1 divergence (5 mismatches) = 36.
        assertEquals(36.0, stats.get("v1").marginMean(), EPSILON);
    }

    private static ViralReference reference()
    {
        List<ViralContig> contigs = List.of(
                new ViralContig("v1", 20, "Virus 1", "Group 1"),
                new ViralContig("v2", 10, "Virus 2", "Group 2"));
        SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("v1", 20),
                new SAMSequenceRecord("v2", 10)));
        return new ViralReference(contigs, dictionary);
    }

    private static SAMFileHeader header()
    {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        header.addSequence(new SAMSequenceRecord("v1", 20));
        header.addSequence(new SAMSequenceRecord("v2", 10));
        return header;
    }

    private static SAMRecord aligned(SAMFileHeader header, String name, int flags, String contig, int start, String cigar, int score)
    {
        return aligned(header, name, flags, contig, start, cigar, score, 0);
    }

    private static SAMRecord aligned(
            SAMFileHeader header, String name, int flags, String contig, int start, String cigar, int score, int editDistance)
    {
        int readLength = TextCigarCodec.decode(cigar).getReadLength();
        SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setFlags(flags);
        record.setReferenceName(contig);
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        record.setReadBases("A".repeat(readLength).getBytes());
        record.setBaseQualities(SamRecordTestUtils.buildDefaultBaseQuals(readLength));
        record.setAttribute(ALIGNMENT_SCORE_ATTRIBUTE, score);
        record.setAttribute(NUM_MUTATONS_ATTRIBUTE, editDistance);
        return record;
    }

    private String writeBam(SAMFileHeader header, List<SAMRecord> records) throws IOException
    {
        File bam = new File(mTempDir.getRoot(), "aligned.bam");
        try(SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, false, bam))
        {
            records.forEach(writer::addAlignment);
        }
        return bam.getPath();
    }
}
