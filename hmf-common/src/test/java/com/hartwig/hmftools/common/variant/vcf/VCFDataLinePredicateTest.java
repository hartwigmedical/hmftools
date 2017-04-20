package com.hartwig.hmftools.common.variant.vcf;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.function.Predicate;

import org.junit.Test;

public class VCFDataLinePredicateTest {

    private final Predicate<String> predicate = new VCFDataLinePredicate();

    @Test
    public void metaInformationLineFails() {
        assertFalse(predicate.test(VCFTestConstants.META_INFORMATION_LINE));
    }

    @Test
    public void headerLineFails() {
        assertFalse(predicate.test(VCFTestConstants.HEADER_LINE));
    }

    @Test
    public void passDataLine1Passes() {
        assertTrue(predicate.test(VCFTestConstants.PASS_DATA_LINE_1));
    }

    @Test
    public void passDataLine2Passes() {
        assertTrue(predicate.test(VCFTestConstants.PASS_DATA_LINE_2));
    }

    @Test
    public void filteredDataLinePasses() {
        assertTrue(predicate.test(VCFTestConstants.FILTERED_DATA_LINE));
    }
}
