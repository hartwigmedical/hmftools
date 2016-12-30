package com.hartwig.hmftools.common.variant.predicate;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.junit.Test;

public class VCFGermlineVariantPredicateTest {

    private static final String NOT_VARIANT = "./.:773,0:773";
    private static final String VARIANT = "1/1:0,3:3:9:103,9,0";

    @Test
    public void checkSNPVariant() {
        final GermlineVariant refSNPVariant = new GermlineVariant(VariantType.SNP, VARIANT, NOT_VARIANT);
        final GermlineVariant noVariant = new GermlineVariant(VariantType.SNP, NOT_VARIANT, NOT_VARIANT);

        VCFGermlineVariantPredicate predicate = new VCFGermlineVariantPredicate(VariantType.SNP, true);
        assertTrue(predicate.test(refSNPVariant));
        assertFalse(predicate.test(noVariant));

        predicate = new VCFGermlineVariantPredicate(VariantType.INDEL, true);
        assertFalse(predicate.test(refSNPVariant));
        assertFalse(predicate.test(noVariant));
    }

    @Test
    public void checkIndelVariant() {
        final GermlineVariant refIndelVariant = new GermlineVariant(VariantType.INDEL, VARIANT, NOT_VARIANT);
        final GermlineVariant noVariant = new GermlineVariant(VariantType.INDEL, NOT_VARIANT, NOT_VARIANT);

        VCFGermlineVariantPredicate predicate = new VCFGermlineVariantPredicate(VariantType.INDEL, true);
        assertTrue(predicate.test(refIndelVariant));
        assertFalse(predicate.test(noVariant));

        predicate = new VCFGermlineVariantPredicate(VariantType.SNP, true);
        assertFalse(predicate.test(refIndelVariant));
        assertFalse(predicate.test(noVariant));
    }
}
