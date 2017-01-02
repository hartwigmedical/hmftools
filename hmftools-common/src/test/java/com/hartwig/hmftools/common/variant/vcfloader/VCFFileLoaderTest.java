package com.hartwig.hmftools.common.variant.vcfloader;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.exception.HealthChecksException;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.junit.Test;

public class VCFFileLoaderTest {

    private static final String VARIANT_PATH = Resources.getResource("variants").getPath();
    private static final String SOMATIC_EXTENSION = "somatics.vcf";

    @Test
    public void canLoadSomaticVCF() throws IOException, HealthChecksException {
        List<SomaticVariant> variants = VCFFileLoader.loadSomaticVCF(VARIANT_PATH, SOMATIC_EXTENSION);
        assertEquals(3, variants.size());
    }
}