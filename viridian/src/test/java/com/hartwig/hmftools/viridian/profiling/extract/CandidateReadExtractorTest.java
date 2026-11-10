package com.hartwig.hmftools.viridian.profiling.extract;

import static java.lang.String.format;
import static java.util.Collections.singleton;

import static com.hartwig.hmftools.viridian.ViridianTestUtils.samFlags;
import static com.hartwig.hmftools.viridian.common.ViridianConstants.VIRAL_READ_EXTRACTION_PARTITION_SIZE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThrows;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

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
import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.cram.CRAIEntry;
import htsjdk.samtools.cram.ref.ReferenceSource;
import htsjdk.samtools.cram.structure.CRAMEncodingStrategy;
import htsjdk.samtools.seekablestream.SeekableFileStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import htsjdk.samtools.util.FileExtensions;

public class CandidateReadExtractorTest
{
    @Rule
    public TemporaryFolder mTempDir = new TemporaryFolder();

    private static final int MIN_SOFT_CLIP_BASES = 20;

    private static final int UNPAIRED = 0;
    private static final int FIRST_OF_PAIR = samFlags(SAMFlag.READ_PAIRED, SAMFlag.FIRST_OF_PAIR);
    private static final int SECOND_OF_PAIR = samFlags(SAMFlag.READ_PAIRED, SAMFlag.SECOND_OF_PAIR);
    private static final int UNMAPPED = samFlags(SAMFlag.READ_UNMAPPED);
    private static final int MATE_UNMAPPED = samFlags(SAMFlag.MATE_UNMAPPED);

    @Test
    public void testExtract() throws IOException
    {
        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, 10000);
        List<SAMRecord> records = List.of(
                placed(header, "plain", UNPAIRED, "chr1", 100, "100M", "AAAAA"),        // no viral signal: dropped
                placed(header, "clip", FIRST_OF_PAIR, "chr1", 150, "20S80M", "CCCCC"),  // soft-clip candidate
                fullyUnmapped(header, "unmap", SECOND_OF_PAIR, "GGGGG"),
                fullyUnmapped(header, "single", UNPAIRED, "TTTTT"));                    // no mate suffix

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
    public void testExtractPartitionBoundary() throws IOException
    {
        // A read overlapping a boundary is sliced by both partitions, but belongs to the one holding its start.

        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, 2 * VIRAL_READ_EXTRACTION_PARTITION_SIZE);
        int start = VIRAL_READ_EXTRACTION_PARTITION_SIZE - 50;   // the 80 aligned bases run past the first partition's end
        List<SAMRecord> records = List.of(placed(header, "spanning", FIRST_OF_PAIR, "chr1", start, "20S80M", "CCCCC"));

        String bam = writeIndexedBam(header, records);
        String fasta = new File(mTempDir.getRoot(), "boundary.fasta").getPath();

        CandidateReadFilter filter = filter(CandidateReadSource.ALL);
        new CandidateReadExtractor(null, filter, 2).extract(bam, fasta);

        assertEquals(List.of(">spanning/1\nCCCCC"), fastaEntries(fasta));
    }

    @Test
    public void testExtractUnmappedReadWithMappedMate() throws IOException
    {
        // Such a read sits in coordinate order, not the fully unmapped block at the end of the BAM. It may be an
        // integration site.

        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, 10000);
        List<SAMRecord> records = List.of(
                placed(header, "unmappedWithMappedMate", FIRST_OF_PAIR | UNMAPPED, "chr1", 500, "*", "TTTTT"));

        String bam = writeIndexedBam(header, records);
        String fasta = new File(mTempDir.getRoot(), "unmapped_with_mapped_mate.fasta").getPath();

        CandidateReadFilter filter = filter(CandidateReadSource.ALL);
        new CandidateReadExtractor(null, filter, 1).extract(bam, fasta);

        assertEquals(List.of(">unmappedWithMappedMate/1\nTTTTT"), fastaEntries(fasta));
    }

    @Test
    public void testExtractFullyUnmappedAndDecoySource() throws IOException
    {
        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, 10000);
        header.addSequence(new SAMSequenceRecord("chrEBV", 10000));
        List<SAMRecord> records = List.of(
                placed(header, "clip", FIRST_OF_PAIR, "chr1", 100, "20S80M", "CCCCC"),
                placed(header, "mappedWithUnmappedMate", FIRST_OF_PAIR | MATE_UNMAPPED, "chr1", 200, "100M", "AAAAA"),
                placed(header, "unmappedWithMappedMate", SECOND_OF_PAIR | UNMAPPED, "chr1", 200, "*", "TTTTT"),
                placed(header, "ebv", FIRST_OF_PAIR, "chrEBV", 300, "100M", "GGGGG"),
                fullyUnmapped(header, "fullyUnmapped", SECOND_OF_PAIR, "ACGTA"));

        String bam = writeIndexedBam(header, records);
        String fasta = new File(mTempDir.getRoot(), "fully_unmapped_and_decoy.fasta").getPath();

        CandidateReadFilter filter = filter(CandidateReadSource.FULLY_UNMAPPED_AND_DECOY);
        new CandidateReadExtractor(null, filter, 2).extract(bam, fasta);

        assertEquals(List.of(">ebv/1\nGGGGG", ">fullyUnmapped/2\nACGTA"), fastaEntries(fasta));
    }

    @Test
    public void testExtractUnindexedBam() throws IOException
    {
        // Sharding the scan by region needs an index.

        SAMFileHeader header = header(SAMFileHeader.SortOrder.unsorted, 10000);
        String bam = writeBam(header, List.of(fullyUnmapped(header, "unmap", UNPAIRED, "GGGGG")));
        String fasta = new File(mTempDir.getRoot(), "candidates.fasta").getPath();

        CandidateReadFilter filter = filter(CandidateReadSource.ALL);
        CandidateReadExtractor extractor = new CandidateReadExtractor(null, filter, 1);

        assertThrows(UserInputError.class, () -> extractor.extract(bam, fasta));
    }

    @Test
    public void testExtractCramUnmappedPartitions() throws IOException
    {
        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, 10000);
        File reference = writeReference(header);
        File cram = writeCram(header, cramRecords(header), reference, 1, null);
        writeCrai(cram);

        assertEquals(List.of(1, 1, 1), slicesPerContainer(cram));
        assertCramCandidates(cram, reference, CRAM_CANDIDATES);
    }

    @Test
    public void testExtractCramMultiSliceContainers() throws IOException
    {
        // Up to two slices per container. The writer gives the mixed slice its own container: [clip, u1],
        // [u2, u3 | u4, u5] and [u6, u7]. With 4 threads each container is its own partition, and the second slice of the
        // middle one must still be read.

        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, 10000);
        List<SAMRecord> records = new ArrayList<>(cramRecords(header));
        records.add(fullyUnmapped(header, "u6", UNPAIRED, "ACGTA"));
        records.add(fullyUnmapped(header, "u7", UNPAIRED, "ACGTA"));
        File reference = writeReference(header);
        File cram = writeCram(header, records, reference, 2, null);
        writeCrai(cram);

        assertEquals(List.of(1, 2, 1), slicesPerContainer(cram));
        List<String> expected = new ArrayList<>(CRAM_CANDIDATES);
        expected.addAll(List.of(">u6\nACGTA", ">u7\nACGTA"));
        assertCramCandidates(cram, reference, expected);
    }

    @Test
    public void testExtractCramMixedContainer() throws IOException
    {
        // [clip, u1] is the only container, so no container holds unmapped reads alone.

        SAMFileHeader header = header(SAMFileHeader.SortOrder.coordinate, 10000);
        List<SAMRecord> records = List.of(
                placed(header, "clip", FIRST_OF_PAIR, "chr1", 100, "20S80M", "C".repeat(100)),
                fullyUnmapped(header, "u1", UNPAIRED, "ACGTA"));
        File reference = writeReference(header);
        File cram = writeCram(header, records, reference, 1, null);
        writeCrai(cram);

        assertEquals(List.of(1), slicesPerContainer(cram));
        assertCramCandidates(cram, reference, List.of(">clip/1\n" + "C".repeat(100), ">u1\nACGTA"));
    }

    @Test
    public void testExtractCramBaiIndex() throws IOException
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
                placed(header, "clip", FIRST_OF_PAIR, "chr1", 100, "20S80M", "C".repeat(100)),
                fullyUnmapped(header, "u1", UNPAIRED, "ACGTA"),
                fullyUnmapped(header, "u2", UNPAIRED, "ACGTA"),
                fullyUnmapped(header, "u3", UNPAIRED, "ACGTA"),
                fullyUnmapped(header, "u4", UNPAIRED, "ACGTA"),
                fullyUnmapped(header, "u5", UNPAIRED, "ACGTA"));
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
        return new CandidateReadFilter(source, MIN_SOFT_CLIP_BASES, singleton("chrEBV"), null);
    }

    private static SAMFileHeader header(SAMFileHeader.SortOrder sortOrder, int contigLength)
    {
        SAMFileHeader header = new SAMFileHeader();
        header.setSortOrder(sortOrder);
        header.addSequence(new SAMSequenceRecord("chr1", contigLength));
        return header;
    }

    // A read with a position: mapped, or unmapped and placed at its mapped mate.
    private static SAMRecord placed(
            SAMFileHeader header, String name, int flags, String contig, int start, String cigar, String bases)
    {
        SAMRecord record = record(header, name, flags, bases);
        record.setReferenceName(contig);
        record.setAlignmentStart(start);
        record.setCigarString(cigar);
        if(!record.getReadUnmappedFlag())
        {
            record.setMappingQuality(60);
        }
        return record;
    }

    // Unmapped, and unpaired or with an unmapped mate, so it sits in the unmapped block at the end of the file.
    private static SAMRecord fullyUnmapped(SAMFileHeader header, String name, int pairFlags, String bases)
    {
        int mateFlags = pairFlags == UNPAIRED ? 0 : MATE_UNMAPPED;
        return record(header, name, pairFlags | UNMAPPED | mateFlags, bases);
    }

    private static SAMRecord record(SAMFileHeader header, String name, int flags, String bases)
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

    // Checks the container layout a test relies on, since the CRAM writer decides it.
    private static List<Integer> slicesPerContainer(File cram) throws IOException
    {
        try(InputStream craiStream = new FileInputStream(cram.getPath() + FileExtensions.CRAM_INDEX))
        {
            Map<Long, Set<Integer>> slicesByContainer = new TreeMap<>();
            for(CRAIEntry entry : CRAMCRAIIndexer.readIndex(craiStream).getCRAIEntries())
            {
                slicesByContainer.computeIfAbsent(entry.getContainerStartByteOffset(), k -> new HashSet<>())
                        .add(entry.getSliceByteOffsetFromCompressionHeaderStart());
            }
            return slicesByContainer.values().stream().map(Set::size).toList();
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
