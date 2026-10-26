package com.hartwig.hmftools.viridian.profiling.extract;

import static java.lang.String.format;
import static java.util.Collections.singleton;

import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_READ_EXTRACTION_PARTITION_SIZE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.viridian.common.UserInputError;

import org.jetbrains.annotations.Nullable;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.CRAMCRAIIndexer;
import htsjdk.samtools.CRAMFileWriter;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.cram.structure.CRAMEncodingStrategy;
import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.FileExtensions;

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
                unmapped(header, "unmap", 0x1 | 0x4 | 0x8 | 0x80, "GGGGG"),           // unmapped candidate, second of pair
                unmapped(header, "single", 0x4, "TTTTT"));                            // unpaired candidate: no mate suffix

        String bam = writeIndexedBam(header, records);

        // Multithreading must not change which reads are candidates, only the output order.
        for(int threads : new int[] { 1, 4 })
        {
            String fasta = new File(mTempDir.getRoot(), "candidates." + threads + ".fasta").getPath();
            CandidateReadFilter filter = filter(CandidateReadSource.ALL);
            new CandidateReadExtractor(null, filter, threads).extract(bam, fasta);

            assertEquals(List.of(">clip/1\nCCCCC", ">single\nTTTTT", ">unmap/2\nGGGGG"), fastaEntries(fasta));
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

        CandidateReadFilter filter = filter(CandidateReadSource.ALL);
        new CandidateReadExtractor(null, filter, 2).extract(bam, fasta);

        assertEquals(List.of(">spanning/1\nCCCCC"), fastaEntries(fasta));
    }

    @Test
    public void testExtractToFastaUnmappedReadWithMappedMateFoundByRegionScan() throws IOException
    {
        // Such a read sits in coordinate order, not the fully unmapped block at the end of the BAM. It may be an integration site.
        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, 10000);
        List<SAMRecord> records = List.of(mapped(header, "unmappedWithMappedMate", 0x1 | 0x4 | 0x40, "chr1", 500, "*", "TTTTT"));

        String bam = writeIndexedBam(header, records);
        String fasta = new File(mTempDir.getRoot(), "unmapped_with_mapped_mate.fasta").getPath();

        CandidateReadFilter filter = filter(CandidateReadSource.ALL);
        new CandidateReadExtractor(null, filter, 1).extract(bam, fasta);

        assertEquals(List.of(">unmappedWithMappedMate/1\nTTTTT"), fastaEntries(fasta));
    }

    @Test
    public void testExtractToFastaFullyUnmappedAndDecoySourceSkipsHostContigs() throws IOException
    {
        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, 10000);
        header.addSequence(new SAMSequenceRecord("chrEBV", 10000));
        List<SAMRecord> records = List.of(
                mapped(header, "clip", 0x1 | 0x40, "chr1", 100, "20S80M", "CCCCC"),
                mapped(header, "mappedWithUnmappedMate", 0x1 | 0x8 | 0x40, "chr1", 200, "100M", "AAAAA"),
                mapped(header, "unmappedWithMappedMate", 0x1 | 0x4 | 0x80, "chr1", 200, "*", "TTTTT"),
                mapped(header, "ebv", 0x1 | 0x40, "chrEBV", 300, "100M", "GGGGG"),
                unmapped(header, "fullyUnmapped", 0x1 | 0x4 | 0x8 | 0x80, "ACGTA"));

        String bam = writeIndexedBam(header, records);
        String fasta = new File(mTempDir.getRoot(), "fully_unmapped_and_decoy.fasta").getPath();

        CandidateReadFilter filter = filter(CandidateReadSource.FULLY_UNMAPPED_AND_DECOY);
        new CandidateReadExtractor(null, filter, 2).extract(bam, fasta);

        assertEquals(List.of(">ebv/1\nGGGGG", ">fullyUnmapped/2\nACGTA"), fastaEntries(fasta));
    }

    @Test
    public void testExtractToFastaUnindexedBamRejected() throws IOException
    {
        // Sharding the scan by region needs an index.
        SAMFileHeader header = header(SAMFileHeader.SortOrder.unsorted, 10000);
        String bam = writeBam(header, List.of(unmapped(header, "unmap", 0x4, "GGGGG")));
        String fasta = new File(mTempDir.getRoot(), "candidates.fasta").getPath();

        CandidateReadFilter filter = filter(CandidateReadSource.ALL);
        CandidateReadExtractor extractor = new CandidateReadExtractor(null, filter, 1);

        assertThrows(UserInputError.class, () -> extractor.extract(bam, fasta));
    }

    @Test
    public void testExtractToFastaCramUnmappedPartitionsWriteEachReadOnce() throws IOException
    {
        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, 10000);
        File reference = writeReference(header);
        File cram = writeCram(header, cramRecords(header), reference, 1, null);
        writeCrai(cram);

        assertCramCandidates(cram, reference, CRAM_CANDIDATES);
    }

    @Test
    public void testExtractToFastaCramMultiSliceContainersReadWhole() throws IOException
    {
        // Up to two slices per container. The writer gives the mixed slice its own container: [clip, u1], [u2, u3 | u4, u5] and
        // [u6, u7]. With 4 threads each container is its own partition, and the second slice of the middle one must still be read.
        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, 10000);
        List<SAMRecord> records = new ArrayList<>(cramRecords(header));
        records.add(unmapped(header, "u6", 0x4, "ACGTA"));
        records.add(unmapped(header, "u7", 0x4, "ACGTA"));
        File reference = writeReference(header);
        File cram = writeCram(header, records, reference, 2, null);
        writeCrai(cram);

        List<String> expected = new ArrayList<>(CRAM_CANDIDATES);
        expected.addAll(List.of(">u6\nACGTA", ">u7\nACGTA"));
        assertCramCandidates(cram, reference, expected);
    }

    @Test
    public void testExtractToFastaCramUnmappedReadOnlyInMixedContainer() throws IOException
    {
        // [clip, u1] is the only container, so no container holds unmapped reads alone.
        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, 10000);
        List<SAMRecord> records = List.of(
                mapped(header, "clip", 0x1 | 0x40, "chr1", 100, "20S80M", "C".repeat(100)),
                unmapped(header, "u1", 0x4, "ACGTA"));
        File reference = writeReference(header);
        File cram = writeCram(header, records, reference, 1, null);
        writeCrai(cram);

        assertCramCandidates(cram, reference, List.of(">clip/1\n" + "C".repeat(100), ">u1\nACGTA"));
    }

    @Test
    public void testExtractToFastaCramWithBaiIndexRejected() throws IOException
    {
        // Only a CRAI lists the unmapped containers, so a CRAM must have one.
        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, 10000);
        File reference = writeReference(header);
        File bai = new File(mTempDir.getRoot(), "reads.cram" + FileExtensions.BAI_INDEX);
        File cram = writeCram(header, cramRecords(header), reference, 1, bai);
        String fasta = new File(mTempDir.getRoot(), "candidates.fasta").getPath();

        CandidateReadFilter filter = filter(CandidateReadSource.ALL);
        CandidateReadExtractor extractor = new CandidateReadExtractor(reference.getPath(), filter, 1);

        assertThrows(UserInputError.class, () -> extractor.extract(cram.getPath(), fasta));
    }

    // Two reads per container: [clip, u1] share a container, then [u2, u3] and [u4, u5].
    // With 4 threads each container is its own unmapped partition. The mapped read sharing the first container must
    // still come only from its region.
    private static List<SAMRecord> cramRecords(SAMFileHeader header)
    {
        return List.of(
                mapped(header, "clip", 0x1 | 0x40, "chr1", 100, "20S80M", "C".repeat(100)),
                unmapped(header, "u1", 0x4, "ACGTA"),
                unmapped(header, "u2", 0x4, "ACGTA"),
                unmapped(header, "u3", 0x4, "ACGTA"),
                unmapped(header, "u4", 0x4, "ACGTA"),
                unmapped(header, "u5", 0x4, "ACGTA"));
    }

    private static final List<String> CRAM_CANDIDATES = List.of(
            ">clip/1\n" + "C".repeat(100), ">u1\nACGTA", ">u2\nACGTA", ">u3\nACGTA", ">u4\nACGTA", ">u5\nACGTA");

    // Each candidate read must appear exactly once, whether the unmapped containers are one partition or several.
    private void assertCramCandidates(File cram, File reference, List<String> expected) throws IOException
    {
        for(int threads : new int[] { 1, 4 })
        {
            String fasta = new File(mTempDir.getRoot(), "cram_candidates." + threads + ".fasta").getPath();
            CandidateReadFilter filter = filter(CandidateReadSource.ALL);
            new CandidateReadExtractor(reference.getPath(), filter, threads).extract(cram.getPath(), fasta);

            assertEquals(expected, fastaEntries(fasta));
        }
    }

    private static CandidateReadFilter filter(CandidateReadSource source)
    {
        return new CandidateReadFilter(source, 20, singleton("chrEBV"), null);
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

    private File writeReference(SAMFileHeader header) throws IOException
    {
        SAMSequenceRecord contig = header.getSequence(0);
        File fasta = new File(mTempDir.getRoot(), "ref.fasta");
        String nameLine = ">" + contig.getSequenceName() + "\n";
        Files.writeString(fasta.toPath(), nameLine + "A".repeat(contig.getSequenceLength()) + "\n");
        Files.writeString(
                new File(fasta.getPath() + FileExtensions.FASTA_INDEX).toPath(),
                format(
                        "%s\t%d\t%d\t%d\t%d\n",
                        contig.getSequenceName(), contig.getSequenceLength(), nameLine.length(),
                        contig.getSequenceLength(), contig.getSequenceLength() + 1));
        return fasta;
    }

    // Two reads per slice, and a mixed-reference slice where fewer than two reads remain for a contig.
    private File writeCram(
            SAMFileHeader header, List<SAMRecord> records, File reference, int slicesPerContainer, @Nullable File baiFile)
            throws IOException
    {
        File cram = new File(mTempDir.getRoot(), "reads.cram");
        CRAMEncodingStrategy encoding = new CRAMEncodingStrategy()
                .setMinimumSingleReferenceSliceSize(2)
                .setReadsPerSlice(2)
                .setSlicesPerContainer(slicesPerContainer);
        try(OutputStream cramStream = new FileOutputStream(cram);
                OutputStream baiStream = baiFile != null ? new FileOutputStream(baiFile) : null)
        {
            CRAMFileWriter writer = new CRAMFileWriter(
                    encoding, cramStream, baiStream, true, new ReferenceSource(reference), header, cram.getName());
            records.forEach(writer::addAlignment);
            writer.close();
        }
        return cram;
    }

    private static void writeCrai(File cram) throws IOException
    {
        try(SeekableStream cramStream = new SeekableFileStream(cram);
                OutputStream craiStream = new FileOutputStream(cram.getPath() + FileExtensions.CRAM_INDEX))
        {
            CRAMCRAIIndexer.writeIndex(cramStream, craiStream);
        }
    }

    private static List<String> fastaEntries(String fasta) throws IOException
    {
        List<String> lines = Files.readAllLines(new File(fasta).toPath());
        List<String> entries = new ArrayList<>();
        for(int i = 0; i < lines.size(); i += 2)
        {
            entries.add(lines.get(i) + "\n" + lines.get(i + 1));
        }
        Collections.sort(entries);
        return entries;
    }
}
