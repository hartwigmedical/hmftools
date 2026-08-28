package com.hartwig.hmftools.virusdetect;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.bwa.IBwaMemAligner;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class ViralReadAlignerTest
{
    @Rule
    public TemporaryFolder mTempDir = new TemporaryFolder();

    private static final String R1_BASES = "A".repeat(100);
    private static final String R2_BASES = "C".repeat(100);
    private static final String R3_BASES = "G".repeat(100);

    // A read is kept on every contig it hits above the min-score fraction (0.5 -> min score 50); per-hit fate inline below.
    @Test
    public void testWritesAlignmentsAboveThresholdWithContigAndScore() throws IOException
    {
        Map<String, List<BwaMemAlignment>> alignments = Map.of(
                R1_BASES, List.of(
                        alignment(0, 0, 9, 80, "100M"),        // forward primary on hpv16, pos 10
                        alignment(0x110, 1, 4, 60, "10H90M"),  // reverse secondary on hpv18, pos 5
                        alignment(0x100, 0, 200, 10, "100M")), // below threshold: dropped
                R2_BASES, List.of(noHit()),                    // no viral alignment: dropped
                R3_BASES, List.of(alignment(0, 0, 0, 40, "100M"))); // only hit below threshold: dropped

        ViralReadAligner aligner = new ViralReadAligner(new FakeAligner(alignments), header(), 0.5, 2);

        String fasta = writeFasta();
        String bam = new File(mTempDir.getRoot(), "aligned.bam").getPath();
        aligner.align(fasta, bam);

        List<SAMRecord> records = readBam(bam);
        assertEquals(2, records.size());

        SAMRecord primary = records.get(0);
        assertEquals("r1", primary.getReadName());
        assertEquals("hpv16", primary.getReferenceName());
        assertEquals(10, primary.getAlignmentStart());
        assertEquals("100M", primary.getCigarString());
        assertEquals(80, (int) primary.getIntegerAttribute("AS"));
        assertFalse(primary.getReadNegativeStrandFlag());
        assertEquals(R1_BASES, primary.getReadString());

        SAMRecord secondary = records.get(1);
        assertEquals("r1", secondary.getReadName());
        assertEquals("hpv18", secondary.getReferenceName());
        assertEquals(5, secondary.getAlignmentStart());
        assertEquals("10H90M", secondary.getCigarString());
        assertEquals(60, (int) secondary.getIntegerAttribute("AS"));
        assertTrue(secondary.getReadNegativeStrandFlag());
        assertEquals("*", secondary.getReadString());
    }

    private static BwaMemAlignment alignment(int samFlag, int refId, int refStart, int score, String cigar)
    {
        return new BwaMemAlignment(samFlag, refId, refStart, refStart + 100, 0, 100, 60, 0, score, 0, cigar, "100", null, -1, -1, 0);
    }

    private static BwaMemAlignment noHit()
    {
        return new BwaMemAlignment(0x4, -1, -1, -1, 0, 0, 0, 0, 0, 0, "*", null, null, -1, -1, 0);
    }

    private static SAMFileHeader header()
    {
        SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("hpv16", 5000),
                new SAMSequenceRecord("hpv18", 5000)));
        SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        return header;
    }

    private String writeFasta() throws IOException
    {
        File fasta = new File(mTempDir.getRoot(), "candidates.fasta");
        String content = ">r1\n" + R1_BASES + "\n>r2\n" + R2_BASES + "\n>r3\n" + R3_BASES + "\n";
        Files.writeString(fasta.toPath(), content);
        return fasta.getPath();
    }

    private static List<SAMRecord> readBam(String bam)
    {
        List<SAMRecord> records = new ArrayList<>();
        try(SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(new File(bam)))
        {
            for(SAMRecord record : reader)
            {
                records.add(record);
            }
        }
        catch(IOException e)
        {
            throw new RuntimeException(e);
        }
        return records;
    }

    private static class FakeAligner implements IBwaMemAligner
    {
        private final Map<String, List<BwaMemAlignment>> mBySequence;

        FakeAligner(Map<String, List<BwaMemAlignment>> bySequence)
        {
            mBySequence = bySequence;
        }

        @Override
        public List<List<BwaMemAlignment>> alignSequences(List<byte[]> sequences)
        {
            return sequences.stream().map(sequence -> mBySequence.get(new String(sequence))).toList();
        }
    }
}
