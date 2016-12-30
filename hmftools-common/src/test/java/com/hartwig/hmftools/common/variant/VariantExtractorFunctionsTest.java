package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class VariantExtractorFunctionsTest {

    @Test
    public void canDetectSNPAndIndel() {
        String ref = "C";
        String singleAltSnp = "G";
        String multipleAltSnp = "G,T";
        String singleIndel = "AT";

        assertEquals(VariantType.SNP, VariantExtractorFunctions.extractVCFType(ref, singleAltSnp));
        assertEquals(VariantType.SNP, VariantExtractorFunctions.extractVCFType(ref, multipleAltSnp));
        assertEquals(VariantType.INDEL, VariantExtractorFunctions.extractVCFType(ref, singleIndel));
    }
}