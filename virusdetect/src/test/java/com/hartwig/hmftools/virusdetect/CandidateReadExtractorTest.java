package com.hartwig.hmftools.virusdetect;

import static java.util.Collections.singleton;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class CandidateReadExtractorTest
{
    @Rule
    public TemporaryFolder mTempDir = new TemporaryFolder();

    // Filtering, dedup, and mate-numbered single-end output exercised together (per-read fate inline below).
    @Test
    public void testWritesFilteredDedupedReadsWithMateSuffix() throws IOException
    {
        SAMFileHeader header = header();
        List<SAMRecord> records = List.of(
                mapped(header, "plain", 0, "chr1", "100M", "AAAAA"),
                unmapped(header, "dup", 0x4 | 0x400, "TTTTT"),                  // duplicate unmapped: dropped
                mapped(header, "clip", 0x1 | 0x40, "chr1", "20S80M", "CCCCC"),  // soft-clip candidate, first of pair
                unmapped(header, "unmap", 0x1 | 0x4 | 0x80, "GGGGG"));          // unmapped candidate, second of pair

        String bam = writeBam(header, records);
        String fasta = new File(mTempDir.getRoot(), "candidates.fasta").getPath();

        int count = new CandidateReadExtractor(null, new CandidateReadFilter(20, singleton("chrEBV"))).extractToFasta(bam, fasta);

        assertEquals(2, count);
        assertEquals(List.of(">clip/1", "CCCCC", ">unmap/2", "GGGGG"), Files.readAllLines(new File(fasta).toPath()));
    }

    private static SAMFileHeader header()
    {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        header.addSequence(new SAMSequenceRecord("chr1", 10000));
        return header;
    }

    private static SAMRecord mapped(SAMFileHeader header, String name, int flags, String contig, String cigar, String bases)
    {
        SAMRecord record = baseRecord(header, name, flags, bases);
        record.setReferenceName(contig);
        record.setAlignmentStart(100);
        record.setCigarString(cigar);
        record.setMappingQuality(60);
        return record;
    }

    private static SAMRecord unmapped(SAMFileHeader header, String name, int flags, String bases)
    {
        return baseRecord(header, name, flags, bases);
    }

    private static SAMRecord baseRecord(SAMFileHeader header, String name, int flags, String bases)
    {
        SAMRecord record = new SAMRecord(header);
        record.setReadName(name);
        record.setFlags(flags);
        record.setReadBases(bases.getBytes());
        record.setBaseQualities(SamRecordTestUtils.buildDefaultBaseQuals(bases.length()));
        return record;
    }

    private String writeBam(SAMFileHeader header, List<SAMRecord> records) throws IOException
    {
        File bam = new File(mTempDir.getRoot(), "reads.bam");
        try(SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, false, bam))
        {
            records.forEach(writer::addAlignment);
        }
        return bam.getPath();
    }
}
