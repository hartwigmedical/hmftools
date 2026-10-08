package com.hartwig.hmftools.viridian.common;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertThrows;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

public class AtomicOutputTest
{
    @Rule
    public TemporaryFolder mTempDir = new TemporaryFolder();

    @Test
    public void testWrite() throws IOException
    {
        String file = path("out.fasta");
        writeText(file, "old");

        String writtenTo = AtomicOutput.write(
                file, tempFile ->
                {
                    writeText(tempFile, "new");
                    return tempFile;
                });

        assertEquals(path("tmp.out.fasta"), writtenTo);
        assertEquals("new", Files.readString(Path.of(file)));
        assertEquals(Set.of("out.fasta"), dirContents());
    }

    @Test
    public void testWriteFailure() throws IOException
    {
        String file = path("out.fasta");
        writeText(file, "old");
        IOException failure = new IOException("writer failed");

        IOException thrown = assertThrows(
                IOException.class, () -> AtomicOutput.write(
                        file, tempFile ->
                        {
                            writeText(tempFile, "partial");
                            throw failure;
                        }));

        assertSame(failure, thrown);
        // The previous output is gone too, so a later reuse cannot pick up a file from an earlier run.
        assertEquals(Set.of(), dirContents());
    }

    @Test
    public void testWriteIndexedBam() throws IOException
    {
        int recordCount = AtomicOutput.writeIndexedBam(path("out.bam"), tempBam -> writeBam(tempBam, false));

        assertEquals(1, recordCount);
        assertEquals(Set.of("out.bam", "out.bam.bai"), dirContents());
    }

    @Test
    public void testWriteIndexedBamFailure() throws IOException
    {
        String bam = path("out.bam");
        AtomicOutput.writeIndexedBam(bam, tempBam -> writeBam(tempBam, false));

        assertThrows(IllegalStateException.class, () -> AtomicOutput.writeIndexedBam(bam, tempBam -> writeBam(tempBam, true)));

        assertEquals(Set.of(), dirContents());
    }

    @Test
    public void testWriteIndexedBamRejectsNonBam()
    {
        assertThrows(IllegalArgumentException.class, () -> AtomicOutput.writeIndexedBam(path("out.sam"), tempBam -> 0));
    }

    private String path(String fileName)
    {
        return new File(mTempDir.getRoot(), fileName).getPath();
    }

    private Set<String> dirContents()
    {
        return Arrays.stream(mTempDir.getRoot().list()).collect(Collectors.toSet());
    }

    private static void writeText(String file, String text) throws IOException
    {
        Files.writeString(Path.of(file), text);
    }

    // Fails after closing, so the index has been written too and both must be cleaned up.
    private static int writeBam(String bamFile, boolean failAfterClose)
    {
        SAMFileHeader header = new SAMFileHeader();
        header.setSequenceDictionary(new SAMSequenceDictionary(List.of(new SAMSequenceRecord("virus", 1000))));
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        try(SAMFileWriter writer = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(header, false, new File(bamFile)))
        {
            SAMRecord record = new SAMRecord(header);
            record.setReadName("read");
            record.setReferenceIndex(0);
            record.setAlignmentStart(1);
            record.setCigarString("4M");
            record.setReadBases("ACGT".getBytes());
            writer.addAlignment(record);
        }

        if(failAfterClose)
        {
            throw new IllegalStateException("writer failed");
        }
        return 1;
    }
}
