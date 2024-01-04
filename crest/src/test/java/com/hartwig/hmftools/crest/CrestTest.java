package com.hartwig.hmftools.crest;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CrestTest
{

    private static final String VCF_FILE = Resources.getResource("purple/sample.purple.germline.vcf.gz").getPath();

    private static final String VCF_DIR = Resources.getResource("purple").getPath();

    private static final String WGS_SAMPLE = "sample";

    private static final String RNA_SAMPLE = "COLO829v003T_RNA";

    private static final String OUTPUT_DIR =
            System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp" + File.separator;

    private static final double EPSILON = 1E-8;

    @Test
    public void shouldFlagMismatchedSample() throws IOException
    {
        String expectedFile = OUTPUT_DIR + File.separator + "sample.CrestCheckFailed";
        cleanCheckFile(expectedFile);

        CrestApplication crestApplication = new CrestApplication(VCF_DIR, OUTPUT_DIR, WGS_SAMPLE, RNA_SAMPLE,
                10, 1, 0.9, false);
        crestApplication.run();
        validateCheckFile(expectedFile);
    }

    @Test
    public void shouldAcceptMatchedSample() throws IOException
    {
        String expectedFile = OUTPUT_DIR + File.separator + "sample.CrestCheckSucceeded";
        cleanCheckFile(expectedFile);

        // rather than cook up another test file, just lower the threshold
        CrestApplication crestApplication = new CrestApplication(VCF_DIR, OUTPUT_DIR, WGS_SAMPLE, RNA_SAMPLE,
                10, 1, 0.3, false);
        crestApplication.run();

        validateCheckFile(expectedFile);
    }

    @Test(expected = RuntimeException.class)
    public void shouldErrorIfNoRnaSample() throws IOException
    {
        CrestApplication crestApplication = new CrestApplication(VCF_DIR, OUTPUT_DIR, WGS_SAMPLE, "WRONG_NAME",
                10, 1, 0.9, true);
        crestApplication.run();
    }

    @Test
    public void shouldComputeCorrectAlleleRatio() throws IOException
    {
        CrestApplication crestApplication = new CrestApplication("purple_dir_not_needed", OUTPUT_DIR,
                "sample_not_needed", RNA_SAMPLE,
                10, 1, 0.9, true);
        double x = crestApplication.computeRnaSupportRatio(VCF_FILE);
        assertEquals(0.47058823, x, EPSILON);
    }

    private void cleanCheckFile(@NotNull String filename)
    {
        File expectedFile = new File(filename);
        if(expectedFile.exists())
        {
            System.out.println("Deleting " + expectedFile.getAbsolutePath());
            expectedFile.delete();
        }
    }

    private void validateCheckFile(@NotNull String filename)
    {
        File expectedFile = new File(filename);
        if(!expectedFile.exists())
        {
            fail("Expected " + expectedFile.getAbsolutePath() + " to exist but file not found");
        }
    }
}
