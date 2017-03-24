package com.hartwig.hmftools.common.variant.vcf;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class VCFMetaInformationLinePredicateTest {

    @Test
    public void vcfMetaInformationLinePredicateWorks() {
        final VCFMetaInformationLinePredicate predicate = new VCFMetaInformationLinePredicate();
        assertFalse(predicate.test(VCFTestConstants.HEADER_LINE));
        assertTrue(predicate.test(VCFTestConstants.META_INFORMATION_LINE));
        assertFalse(predicate.test(VCFTestConstants.PASS_DATA_LINE_1));
    }
}