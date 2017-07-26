package com.hartwig.hmftools.common.variant.vcf;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.junit.Test;

public class VCFFileLoaderTest {

    private static final String VARIANT_PATH = Resources.getResource("variants").getPath();
    private static final String SOMATIC_EXTENSION = "somatics.vcf";
    private static final String SOMATIC_GERMLINE_EXTENSION = "somatic_germline.vcf";
    private static final String SINGLE_SAMPLE_EXTENSION = "single_sample_germline.vcf";

    @Test
    public void canLoadSomaticVCFFromBasePath() throws IOException, HartwigException {
        final VCFSomaticFile variantFile = VCFFileLoader.loadSomaticVCF(VARIANT_PATH, SOMATIC_EXTENSION);
        assertEquals("sample", variantFile.sample());
        assertEquals(3, variantFile.variants().size());
    }

    @Test
    public void canLoadSomaticVCFFromFile() throws IOException, HartwigException {
        final String file = VARIANT_PATH + File.separator + SOMATIC_EXTENSION;
        final VCFSomaticFile variantFile = VCFFileLoader.loadSomaticVCF(file);
        assertEquals("sample", variantFile.sample());
        assertEquals(3, variantFile.variants().size());
    }

    @Test
    public void canLoadSomaticGermlineVCF() throws IOException, HartwigException {
        final VCFGermlineFile variantFile = VCFFileLoader.loadGermlineVCF(VARIANT_PATH, SOMATIC_GERMLINE_EXTENSION);
        assertEquals("sampleR", variantFile.refSample());
        assertEquals("sampleT", variantFile.tumorSample());
        assertEquals(3, variantFile.variants().size());
    }

    @Test
    public void canLoadSingleSampleGermlineVCF() throws IOException, HartwigException {
        final VCFGermlineFile variantFile = VCFFileLoader.loadGermlineVCF(VARIANT_PATH, SINGLE_SAMPLE_EXTENSION);
        assertEquals("sample", variantFile.refSample());
        assertNull(variantFile.tumorSample());
        assertEquals(3, variantFile.variants().size());
    }
}