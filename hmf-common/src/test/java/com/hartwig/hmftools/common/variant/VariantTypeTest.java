package com.hartwig.hmftools.common.variant;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class VariantTypeTest {

    @Test
    public void canDetermineSNPAndIndel() {
        final String ref = "C";
        final String singleAltSnp = "G";
        final String multipleAltSnp = "G,T";
        final String singleIndel = "AT";
        final String snpAndIndel = "G,AG";

        assertEquals(VariantType.SNP, VariantType.fromRefAlt(ref, singleAltSnp));
        assertEquals(VariantType.SNP, VariantType.fromRefAlt(ref, multipleAltSnp));
        assertEquals(VariantType.INDEL, VariantType.fromRefAlt(ref, singleIndel));
        assertEquals(VariantType.INDEL, VariantType.fromRefAlt(ref, snpAndIndel));
    }
}