package com.hartwig.hmftools.common.variant.vcfloader;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.function.Predicate;

import org.junit.Test;

public class VCFDataLinePredicateTest {

    @Test
    public void commentLineFails() {
        Predicate<String> predicate = new VCFDataLinePredicate();
        assertFalse(predicate.test(VCFTestConstants.COMMENT_LINE));
    }

    @Test
    public void headerLineFails() {
        Predicate<String> predicate = new VCFDataLinePredicate();
        assertFalse(predicate.test(VCFTestConstants.HEADER_LINE));
    }

    @Test
    public void passDataLine1Passes() {
        Predicate<String> predicate = new VCFDataLinePredicate();
        assertTrue(predicate.test(VCFTestConstants.PASS_DATA_LINE_1));
    }

    @Test
    public void passDataLine2Passes() {
        VCFDataLinePredicate predicate = new VCFDataLinePredicate();
        assertTrue(predicate.test(VCFTestConstants.PASS_DATA_LINE_2));
    }

    @Test
    public void filteredDataLinePasses() {
        Predicate<String> predicate = new VCFDataLinePredicate();
        assertTrue(predicate.test(VCFTestConstants.FILTERED_DATA_LINE));
    }
}
