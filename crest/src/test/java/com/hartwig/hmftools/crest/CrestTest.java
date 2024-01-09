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
    private static final String VCF_FILE = Resources.getResource("purple/tumor_sample.purple.germline.vcf.gz").getPath();
    private static final String VCF_DIR = Resources.getResource("purple").getPath();
    private static final String WGS_SAMPLE = "tumor_sample";
    private static final String RNA_SAMPLE = "rna_sample";

    private static final String OUTPUT_DIR =
            System.getProperty("user.home") + File.separator + "hmf" + File.separator + "tmp" + File.separator;

    private static final double EPSILON = 1E-8;

    @Test
    public void shouldFlagMismatchedSample() throws IOException
    {
        String expectedFile = OUTPUT_DIR + File.separator + "tumor_sample.CrestCheckFailed";
        cleanCheckFile(expectedFile);

        CrestAlgo crestAlgo = new CrestAlgo(VCF_DIR, OUTPUT_DIR, WGS_SAMPLE, RNA_SAMPLE,
                10, 1, 0.9, false);
        crestAlgo.run();
        validateCheckFile(expectedFile);
    }

    @Test
    public void shouldAcceptMatchedSample() throws IOException
    {
        String expectedFile = OUTPUT_DIR + File.separator + "tumor_sample.CrestCheckSucceeded";
        cleanCheckFile(expectedFile);

        // rather than cook up another test file, just lower the threshold
        CrestAlgo crestAlgo = new CrestAlgo(VCF_DIR, OUTPUT_DIR, WGS_SAMPLE, RNA_SAMPLE,
                10, 1, 0.3, false);
        crestAlgo.run();

        validateCheckFile(expectedFile);
    }

    @Test(expected = RuntimeException.class)
    public void shouldErrorIfNoRnaSample() throws IOException
    {
        CrestAlgo crestAlgo = new CrestAlgo(VCF_DIR, OUTPUT_DIR, WGS_SAMPLE, "WRONG_NAME",
                10, 1, 0.9, true);
        crestAlgo.run();
    }

    @Test
    public void shouldComputeCorrectAlleleRatio() throws IOException
    {
        CrestAlgo crestAlgo = new CrestAlgo("purple_dir_not_needed", null, WGS_SAMPLE, RNA_SAMPLE,
                10, 1, 0.9, true);
        double x = crestAlgo.computeRnaSupportRatio(VCF_FILE);
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
