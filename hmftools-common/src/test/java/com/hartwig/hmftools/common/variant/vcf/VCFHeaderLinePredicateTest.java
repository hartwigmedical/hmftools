package com.hartwig.hmftools.common.variant.vcf;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class VCFHeaderLinePredicateTest {

    @Test
    public void vcfHeaderLinePredicateWorks() {
        final VCFHeaderLinePredicate predicate = new VCFHeaderLinePredicate();
        assertTrue(predicate.test(VCFTestConstants.HEADER_LINE));
        assertFalse(predicate.test(VCFTestConstants.META_INFORMATION_LINE));
        assertFalse(predicate.test(VCFTestConstants.PASS_DATA_LINE_1));
    }
}