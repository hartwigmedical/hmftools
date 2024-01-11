package com.hartwig.hmftools.crest;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Test;

public class CrestTest
{
    private static final String MINIMAL_VCF_FILE = Resources.getResource("minimal.vcf").getPath();
    private static final String PURPLE_DIR = Resources.getResource("purple").getPath();
    private static final String WGS_SAMPLE = "tumor_sample";
    private static final String RNA_SAMPLE = "rna_sample";

    private static final double EPSILON = 1E-8;

    @Test
    public void shouldFlagMismatchedSample() throws IOException
    {
        CrestAlgo crestAlgo = new CrestAlgo(PURPLE_DIR, null, WGS_SAMPLE, RNA_SAMPLE,
                10, 1, 0.9, true);
        assertFalse(crestAlgo.crestCheck(MINIMAL_VCF_FILE));
    }

    @Test
    public void shouldAcceptMatchedSample() throws IOException
    {
        // rather than cook up another test file, just lower the threshold
        CrestAlgo crestAlgo = new CrestAlgo(PURPLE_DIR, null, WGS_SAMPLE, RNA_SAMPLE,
                10, 1, 0.3, true);
        assertTrue(crestAlgo.crestCheck(MINIMAL_VCF_FILE));
    }

    @Test(expected = RuntimeException.class)
    public void shouldErrorIfNoRnaSample() throws IOException
    {
        CrestAlgo crestAlgo = new CrestAlgo(PURPLE_DIR, null, WGS_SAMPLE, "WRONG_NAME",
                10, 1, 0.9, true);
        crestAlgo.crestCheck(MINIMAL_VCF_FILE);
    }

    @Test
    public void shouldRunOnRealVcf() throws IOException
    {
        CrestAlgo crestAlgo = new CrestAlgo(PURPLE_DIR, null, WGS_SAMPLE, RNA_SAMPLE,
                10, 1, 0.9, true);
        crestAlgo.run();
    }

    @Test
    public void shouldIgnoreInvalidRecords() throws IOException
    {
        CrestAlgo crestAlgo = new CrestAlgo(PURPLE_DIR, null, WGS_SAMPLE, RNA_SAMPLE,
                10, 1, 0.9, true);
        double x = crestAlgo.computeRnaSupportRatio(MINIMAL_VCF_FILE);
        assertEquals(0.5, x, EPSILON);
    }
}
