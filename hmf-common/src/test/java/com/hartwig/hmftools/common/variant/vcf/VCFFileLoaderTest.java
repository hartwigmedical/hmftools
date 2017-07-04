package com.hartwig.hmftools.common.variant.vcf;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;
import org.junit.Test;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.assertEquals;

public class VCFFileLoaderTest {

    private static final String VARIANT_PATH = Resources.getResource("variants").getPath();
    private static final String SOMATIC_EXTENSION = "somatics.vcf";
    private static final String GERMLINE_EXTENSION = "germline.vcf";

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
    public void canLoadGermlineVCF() throws IOException, HartwigException {
        final VCFGermlineFile variantFile = VCFFileLoader.loadGermlineVCF(VARIANT_PATH, GERMLINE_EXTENSION);
        assertEquals("sampleR", variantFile.refSample());
        assertEquals("sampleT", variantFile.tumorSample());
        assertEquals(3, variantFile.variants().size());
    }
    
    @Test
    public void canLoadGermlineVCFFromFile() throws IOException, HartwigException {
        final String file = VARIANT_PATH + File.separator + GERMLINE_EXTENSION;
        final VCFGermlineFile variantFile = VCFFileLoader.loadGermlineVCF(file);
        assertEquals("sampleR", variantFile.refSample());
        assertEquals("sampleT", variantFile.tumorSample());
        assertEquals(3, variantFile.variants().size());
    }
}