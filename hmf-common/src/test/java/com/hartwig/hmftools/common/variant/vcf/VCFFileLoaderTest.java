package com.hartwig.hmftools.common.variant.vcf;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.io.IOException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;

import org.junit.Test;

public class VCFFileLoaderTest {

    private static final String VARIANT_PATH = Resources.getResource("variants").getPath();
    private static final String SOMATIC_EXTENSION = "somatics.vcf";

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
}