package com.hartwig.hmftools.viridian.detection.read_extract;

import static java.util.Collections.singleton;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_READ_EXTRACTION_PARTITION_SIZE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.viridian.common.UserInputError;

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

    @Test
    public void testExtractToFastaWritesFilteredCandidateReadsWithMateSuffix() throws IOException
    {
        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, 10000);
        List<SAMRecord> records = List.of(
                mapped(header, "plain", 0, "chr1", 100, "100M", "AAAAA"),             // no viral signal: dropped
                mapped(header, "clip", 0x1 | 0x40, "chr1", 150, "20S80M", "CCCCC"),   // soft-clip candidate, first of pair
                unmapped(header, "dup", 0x4 | 0x400, "TTTTT"),                        // duplicate unmapped: dropped
                unmapped(header, "unmap", 0x1 | 0x4 | 0x80, "GGGGG"),                 // unmapped candidate, second of pair
                unmapped(header, "single", 0x4, "TTTTT"));                            // unpaired candidate: no mate suffix

        String bam = writeIndexedBam(header, records);

        // Multithreading must not change which reads are candidates, only the output order.
        for(int threads : new int[] { 1, 4 })
        {
            String fasta = new File(mTempDir.getRoot(), "candidates." + threads + ".fasta").getPath();
            int count = new CandidateReadExtractor(null, new CandidateReadFilter(20, singleton("chrEBV")), threads)
                    .extractToFasta(bam, fasta);

            assertEquals(3, count);
            assertEquals(Set.of(">clip/1\nCCCCC", ">unmap/2\nGGGGG", ">single\nTTTTT"), fastaEntries(fasta));
        }
    }

    @Test
    public void testExtractToFastaReadSpanningPartitionBoundaryWrittenOnce() throws IOException
    {
        // A read overlapping a boundary is sliced by both partitions, but belongs to the one holding its start.
        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, 2 * VIRAL_READ_EXTRACTION_PARTITION_SIZE);
        int start = VIRAL_READ_EXTRACTION_PARTITION_SIZE - 50;   // the 80 aligned bases run past the first partition's end
        List<SAMRecord> records = List.of(mapped(header, "spanning", 0x1 | 0x40, "chr1", start, "20S80M", "CCCCC"));

        String bam = writeIndexedBam(header, records);
        String fasta = new File(mTempDir.getRoot(), "boundary.fasta").getPath();

        int count = new CandidateReadExtractor(null, new CandidateReadFilter(20, singleton("chrEBV")), 2)
                .extractToFasta(bam, fasta);

        assertEquals(1, count);
        assertEquals(Set.of(">spanning/1\nCCCCC"), fastaEntries(fasta));
    }

    @Test
    public void testExtractToFastaPlacedUnmappedReadFoundByRegionScan() throws IOException
    {
        // Such a read sits in coordinate order, not the unplaced block at the end of the BAM. It may be an integration site.
        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, 10000);
        List<SAMRecord> records = List.of(mapped(header, "placed", 0x1 | 0x4 | 0x40, "chr1", 500, "*", "TTTTT"));

        String bam = writeIndexedBam(header, records);
        String fasta = new File(mTempDir.getRoot(), "placed.fasta").getPath();

        int count = new CandidateReadExtractor(null, new CandidateReadFilter(20, singleton("chrEBV")))
                .extractToFasta(bam, fasta);

        assertEquals(1, count);
        assertEquals(Set.of(">placed/1\nTTTTT"), fastaEntries(fasta));
    }

    @Test
    public void testExtractToFastaUnindexedBamRejected() throws IOException
    {
        // Sharding the scan by region needs an index.
        SAMFileHeader header = header(SAMFileHeader.SortOrder.unsorted, 10000);
        String bam = writeBam(header, List.of(unmapped(header, "unmap", 0x4, "GGGGG")));
        String fasta = new File(mTempDir.getRoot(), "candidates.fasta").getPath();

        CandidateReadExtractor extractor = new CandidateReadExtractor(null, new CandidateReadFilter(20, singleton("chrEBV")));

        assertThrows(UserInputError.class, () -> extractor.extractToFasta(bam, fasta));
    }

    private static SAMFileHeader header(SAMFileHeader.SortOrder sortOrder, int contigLength)
    {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(sortOrder);
        header.addSequence(new SAMSequenceRecord("chr1", contigLength));
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
