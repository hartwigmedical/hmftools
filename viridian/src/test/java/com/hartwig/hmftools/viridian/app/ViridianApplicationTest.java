package com.hartwig.hmftools.viridian.app;

import static com.hartwig.hmftools.viridian.app.ViridianApplication.checkVirusDecoyContigsInTumorBam;

import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.io.File;
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
    public void testCheckVirusDecoyContigsInTumorBamAcceptsContigsInHeader()
    {
        String bam = writeHeaderOnlyBam("chr1", "chrEBV", "chrHPV16");

        checkVirusDecoyContigsInTumorBam(bam, null, Set.of());
        checkVirusDecoyContigsInTumorBam(bam, null, Set.of("chrHPV16"));
    }

    @Test
    public void testCheckVirusDecoyContigsInTumorBamRejectsMissingContigs()
    {
        String bam = writeHeaderOnlyBam("chr1", "chrEBV");

        UserInputError error = assertThrows(
                UserInputError.class, () -> checkVirusDecoyContigsInTumorBam(bam, null, Set.of("chrEBV", "chrHPV16", "chrHBV")));

        assertTrue(error.getMessage().endsWith(": chrHBV, chrHPV16"));
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
