package com.hartwig.hmftools.common.variant.predicate;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.variant.VCFGermlineData;
import com.hartwig.hmftools.common.variant.VCFType;
import com.hartwig.hmftools.common.variant.predicate.VCFGermlineVariantPredicate;

import org.junit.Test;

public class VCFGermlineVariantPredicateTest {

    private static final String NOT_VARIANT = "./.:773,0:773";
    private static final String VARIANT = "1/1:0,3:3:9:103,9,0";

    @Test
    public void checkSNPVariant() {
        final VCFGermlineData vcfGermlineVariantData = new VCFGermlineData(VCFType.SNP, VARIANT, NOT_VARIANT);
        final VCFGermlineData vcfGermlineNotVariantData = new VCFGermlineData(VCFType.SNP, NOT_VARIANT, NOT_VARIANT);

        VCFGermlineVariantPredicate predicate = new VCFGermlineVariantPredicate(VCFType.SNP, true);
        assertTrue(predicate.test(vcfGermlineVariantData));
        assertFalse(predicate.test(vcfGermlineNotVariantData));

        predicate = new VCFGermlineVariantPredicate(VCFType.INDELS, true);
        assertFalse(predicate.test(vcfGermlineVariantData));
        assertFalse(predicate.test(vcfGermlineNotVariantData));
    }

    @Test
    public void checkIndelVariant() {
        final VCFGermlineData vcfGermlineVariantData = new VCFGermlineData(VCFType.INDELS, VARIANT, NOT_VARIANT);
        final VCFGermlineData vcfGermlineNotVariantData = new VCFGermlineData(VCFType.INDELS, NOT_VARIANT, NOT_VARIANT);

        VCFGermlineVariantPredicate predicate = new VCFGermlineVariantPredicate(VCFType.INDELS, true);
        assertTrue(predicate.test(vcfGermlineVariantData));
        assertFalse(predicate.test(vcfGermlineNotVariantData));

        predicate = new VCFGermlineVariantPredicate(VCFType.SNP, true);
        assertFalse(predicate.test(vcfGermlineVariantData));
        assertFalse(predicate.test(vcfGermlineNotVariantData));
    }
}
