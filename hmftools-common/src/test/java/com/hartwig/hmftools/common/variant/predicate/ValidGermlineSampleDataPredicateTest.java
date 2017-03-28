package com.hartwig.hmftools.common.variant.predicate;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.variant.GermlineSampleData;

import org.junit.Test;

public class ValidGermlineSampleDataPredicateTest {

    private static final GermlineSampleData VALID = new GermlineSampleData("1/1", 10, 5);
    private static final GermlineSampleData INVALID = new GermlineSampleData("./.", 10, 5);

    @Test
    public void canCheckGermlineSampleDataCorrectly() {
        final Predicate<GermlineSampleData> predicate = new ValidGermlineSampleDataPredicate();
        assertTrue(predicate.test(VALID));
        assertFalse(predicate.test(INVALID));
    }
}
