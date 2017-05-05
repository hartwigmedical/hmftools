package com.hartwig.hmftools.common.variant.predicate;

import com.hartwig.hmftools.common.variant.GermlineSampleData;
import com.hartwig.hmftools.common.variant.ImmutableGermlineSampleData;
import org.junit.Test;

import java.util.function.Predicate;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class ValidGermlineSampleDataPredicateTest {

    private static final GermlineSampleData VALID = ImmutableGermlineSampleData.of("1/1", 10, 5, 15);
    private static final GermlineSampleData INVALID = ImmutableGermlineSampleData.of("./.", 10, 5, 15);

    @Test
    public void canCheckGermlineSampleDataCorrectly() {
        final Predicate<GermlineSampleData> predicate = new ValidGermlineSampleDataPredicate();
        assertTrue(predicate.test(VALID));
        assertFalse(predicate.test(INVALID));
    }
}
