package com.hartwig.hmftools.common.variant.vcfloader;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.junit.Test;

public class VCFFileLoaderTest {

    private static final String VARIANT_PATH = Resources.getResource("variants").getPath();
    private static final String SOMATIC_EXTENSION = "somatics.vcf";
    private static final String GERMLINE_EXTENSION = "germline.vcf";

    @Test
    public void canLoadSomaticVCF() throws IOException, HartwigException {
        List<SomaticVariant> variants = VCFFileLoader.loadSomaticVCF(VARIANT_PATH, SOMATIC_EXTENSION);
        assertEquals(3, variants.size());
    }

    @Test
    public void canLoadGermlineVCF() throws IOException, HartwigException {
        List<GermlineVariant> variants = VCFFileLoader.loadGermlineVCF(VARIANT_PATH, GERMLINE_EXTENSION);
        assertEquals(3, variants.size());
    }
}