package com.hartwig.hmftools.virusdetect;

import static java.util.Collections.singleton;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

public class CandidateReadExtractorTest
{
    @Rule
    public TemporaryFolder mTempDir = new TemporaryFolder();

    // Filtering, dedup, and mate-numbered single-end output exercised together (per-read fate inline below).
    // Run single- and multi-threaded: sharding must not change which reads are candidates, only the output order.
    @Test
    public void testWritesFilteredDedupedReadsWithMateSuffix() throws IOException
    {
        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate);
        List<SAMRecord> records = List.of(
                mapped(header, "plain", 0, "chr1", 100, "100M", "AAAAA"),             // no viral signal: dropped
                mapped(header, "clip", 0x1 | 0x40, "chr1", 150, "20S80M", "CCCCC"),   // soft-clip candidate, first of pair
                unmapped(header, "dup", 0x4 | 0x400, "TTTTT"),                        // duplicate unmapped: dropped
                unmapped(header, "unmap", 0x1 | 0x4 | 0x80, "GGGGG"));                // unmapped candidate, second of pair

        String bam = writeIndexedBam(header, records);

        for(int threads : new int[] { 1, 4 })
        {
            String fasta = new File(mTempDir.getRoot(), "candidates." + threads + ".fasta").getPath();
            int count = new CandidateReadExtractor(null, new CandidateReadFilter(20, singleton("chrEBV")), threads)
                    .extractToFasta(bam, fasta);

            assertEquals(2, count);
            assertEquals(Set.of(">clip/1\nCCCCC", ">unmap/2\nGGGGG"), fastaEntries(fasta));
        }
    }

    // Sharding the scan by region needs an index, so an unindexed BAM is rejected rather than silently handled.
    @Test
    public void testUnindexedBamRejected() throws IOException
    {
        SAMFileHeader header = header(SAMFileHeader.SortOrder.unsorted);
        String bam = writeBam(header, List.of(unmapped(header, "unmap", 0x4, "GGGGG")));
        String fasta = new File(mTempDir.getRoot(), "candidates.fasta").getPath();

        CandidateReadExtractor extractor = new CandidateReadExtractor(null, new CandidateReadFilter(20, singleton("chrEBV")));

        assertThrows(UserInputError.class, () -> extractor.extractToFasta(bam, fasta));
    }

    private static SAMFileHeader header(SAMFileHeader.SortOrder sortOrder)
    {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(sortOrder);
        header.addSequence(new SAMSequenceRecord("chr1", 10000));
        return header;
    }

    private static SAMRecord mapped(SAMFileHeader header, String name, int flags, String contig, int start, String cigar, String bases)
    {
        SAMRecord record = baseRecord(header, name, flags, bases);
        record.setReferenceName(contig);
        record.setAlignmentStart(start);
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

    private String writeIndexedBam(SAMFileHeader header, List<SAMRecord> records) throws IOException
    {
        File bam = new File(mTempDir.getRoot(), "reads.sorted.bam");
        try(SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, true, bam))
        {
            records.forEach(writer::addAlignment);
        }
        return bam.getPath();
    }

    private static Set<String> fastaEntries(String fasta) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(fasta).toPath());
        Set<String> entries = new HashSet<>();
        for(int i = 0; i < lines.size(); i += 2)
        {
            entries.add(lines.get(i) + "\n" + lines.get(i + 1));
        }
        return entries;
    }
}
