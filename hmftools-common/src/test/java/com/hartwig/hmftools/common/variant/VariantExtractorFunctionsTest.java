package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class VariantExtractorFunctionsTest {

    @Test
    public void canDetermineSNPAndIndel() {
        final String ref = "C";
        final String singleAltSnp = "G";
        final String multipleAltSnp = "G,T";
        final String singleIndel = "AT";
        final String snpAndIndel = "G,AG";

        assertEquals(VariantType.SNP, VariantExtractorFunctions.determineVariantType(ref, singleAltSnp));
        assertEquals(VariantType.SNP, VariantExtractorFunctions.determineVariantType(ref, multipleAltSnp));
        assertEquals(VariantType.INDEL, VariantExtractorFunctions.determineVariantType(ref, singleIndel));
        assertEquals(VariantType.INDEL, VariantExtractorFunctions.determineVariantType(ref, snpAndIndel));
    }
}