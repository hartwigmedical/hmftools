package com.hartwig.hmftools.viridian.app;

import static com.hartwig.hmftools.viridian.app.ViridianApplication.checkVirusDecoyContigsInTumorBam;

import static org.junit.Assert.assertThrows;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Set;

import com.hartwig.hmftools.viridian.common.UserInputError;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMSequenceRecord;

public class ViridianApplicationTest
{
    @Rule
    public TemporaryFolder mTempDir = new TemporaryFolder();

    @Test
    public void testCheckVirusDecoyContigsInTumorBam() throws IOException
    {
        String bam = writeHeaderOnlyBam("chr1", "chrEBV", "chrHPV16");

        String reference = writeReference();

        checkVirusDecoyContigsInTumorBam(bam, reference, Set.of());
        checkVirusDecoyContigsInTumorBam(bam, reference, Set.of("chrHPV16"));
    }

    @Test
    public void testCheckVirusDecoyContigsInTumorBamMissingContigs() throws IOException
    {
        String bam = writeHeaderOnlyBam("chr1", "chrEBV");
        String reference = writeReference();

        assertThrows(
                UserInputError.class,
                () -> checkVirusDecoyContigsInTumorBam(bam, reference, Set.of("chrEBV", "chrHPV16", "chrHBV")));
    }

    // Read only for CRAM, but the reader opens it for BAM too.
    private String writeReference() throws IOException
    {
        File fasta = new File(mTempDir.getRoot(), "ref.fasta");
        Files.writeString(fasta.toPath(), ">chr1\nACGT\n");
        return fasta.getPath();
    }

    private String writeHeaderOnlyBam(String... contigs)
    {
        SAMFileHeader header = new SAMFileHeader();
        for(String contig : contigs)
        {
            header.addSequence(new SAMSequenceRecord(contig, 10000));
        }

        File bam = new File(mTempDir.getRoot(), "reads.bam");
        new SAMFileWriterFactory().makeBAMWriter(header, false, bam).close();
        return bam.getPath();
    }
}
