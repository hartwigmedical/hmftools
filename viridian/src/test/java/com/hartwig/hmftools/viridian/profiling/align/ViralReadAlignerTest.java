package com.hartwig.hmftools.viridian.profiling.align;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.ALIGNMENT_SCORE_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.MISMATCHES_AND_DELETIONS_ATTRIBUTE;
import static com.hartwig.hmftools.common.bam.SamRecordUtils.NUM_MUTATONS_ATTRIBUTE;
import static com.hartwig.hmftools.viridian.ViridianTestUtils.noHit;
import static com.hartwig.hmftools.viridian.ViridianTestUtils.samFlags;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.bwa.MockBwaMemAligner;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.ValidationStringency;

public class ViralReadAlignerTest
{
    @Rule
    public TemporaryFolder mTempDir = new TemporaryFolder();

    // Three reads, so alignment spans two chunks.
    private static final int CHUNK_SIZE = 2;

    private static final int HPV16_REF_ID = 0;
    private static final int HPV18_REF_ID = 1;

    private static final String R1_BASES = "A".repeat(100);
    private static final String R2_BASES = "C".repeat(100);
    private static final String R3_BASES = "A".repeat(50) + "C".repeat(50);

    @Test
    public void testAlignWritesEveryViralAlignmentWithContigAndScore() throws IOException
    {
        // The aligner's own minimum-score floor sits upstream, so low-scoring alignments must still be written here.

        int secondary = samFlags(SAMFlag.SECONDARY_ALIGNMENT);
        int reverseSecondary = samFlags(SAMFlag.SECONDARY_ALIGNMENT, SAMFlag.READ_REVERSE_STRAND);
        Map<String, List<BwaMemAlignment>> alignments = Map.of(
                R1_BASES, List.of(
                        alignment(0, HPV16_REF_ID, 9, 80, "100M"),
                        alignment(reverseSecondary, HPV18_REF_ID, 4, 60, "10H90M"),
                        alignment(secondary, HPV16_REF_ID, 200, 10, "100M")),
                R2_BASES, List.of(noHit()),
                R3_BASES, List.of(
                        alignment(0, HPV16_REF_ID, 0, 40, "100M"),
                        alignment(reverseSecondary, HPV18_REF_ID, 49, 40, "100M")));

        ViralReadAligner aligner = new ViralReadAligner(MockBwaMemAligner.forQueries(alignments), header(), CHUNK_SIZE);

        String fasta = writeFasta();
        String bam = new File(mTempDir.getRoot(), "aligned.bam").getPath();
        aligner.align(fasta, bam);

        // Coordinate sorted, so ordered by contig then position. A hard-clipped alignment omits its bases. A reverse
        // strand alignment holds the reverse complement of the read's bases.
        List<WrittenRecord> expected = List.of(
                new WrittenRecord("r3", 0, "hpv16", 1, "100M", R3_BASES, 40, 0, "100"),
                new WrittenRecord("r1", 0, "hpv16", 10, "100M", R1_BASES, 80, 0, "100"),
                new WrittenRecord("r1", secondary, "hpv16", 201, "100M", R1_BASES, 10, 0, "100"),
                new WrittenRecord("r1", reverseSecondary, "hpv18", 5, "10H90M", "*", 60, 0, "90"),
                new WrittenRecord("r3", reverseSecondary, "hpv18", 50, "100M", "G".repeat(50) + "T".repeat(50), 40, 0, "100"));
        assertEquals(expected, readBam(bam));
    }

    private record WrittenRecord(
            String readName, int flags, String contig, int alignmentStart, String cigar, String bases, int alignerScore,
            int editDistance, String mismatchesAndDeletions)
    {
        static WrittenRecord from(SAMRecord record)
        {
            return new WrittenRecord(
                    record.getReadName(), record.getFlags(), record.getReferenceName(), record.getAlignmentStart(),
                    record.getCigarString(), record.getReadString(), record.getIntegerAttribute(ALIGNMENT_SCORE_ATTRIBUTE),
                    record.getIntegerAttribute(NUM_MUTATONS_ATTRIBUTE), record.getStringAttribute(MISMATCHES_AND_DELETIONS_ATTRIBUTE));
        }
    }

    // BWA reference positions are 0-based. Every alignment has no edits, so NM is 0 and MD is the aligned length.
    private static BwaMemAlignment alignment(int samFlag, int refId, int refStart, int score, String cigar)
    {
        String mdTag = String.valueOf(TextCigarCodec.decode(cigar).getReferenceLength());
        return new BwaMemAlignment(
                samFlag, refId, refStart, refStart + 100, 0, 100, 60, 0, score, 0, cigar, mdTag, null, -1, -1, 0);
    }

    private static SAMFileHeader header()
    {
        SAMSequenceDictionary dictionary = new SAMSequenceDictionary(List.of(
                new SAMSequenceRecord("hpv16", 5000),
                new SAMSequenceRecord("hpv18", 5000)));
        SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(dictionary);
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
        return header;
    }

    private String writeFasta() throws IOException
    {
        File fasta = new File(mTempDir.getRoot(), "candidates.fasta");
        String content = ">r1\n" + R1_BASES + "\n>r2\n" + R2_BASES + "\n>r3\n" + R3_BASES + "\n";
        Files.writeString(fasta.toPath(), content);
        return fasta.getPath();
    }

    private static List<WrittenRecord> readBam(String bam) throws IOException
    {
        List<WrittenRecord> records = new ArrayList<>();
        try(SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.STRICT).open(new File(bam)))
        {
            for(SAMRecord record : reader)
            {
                records.add(WrittenRecord.from(record));
            }
        }
        return records;
    }
}
