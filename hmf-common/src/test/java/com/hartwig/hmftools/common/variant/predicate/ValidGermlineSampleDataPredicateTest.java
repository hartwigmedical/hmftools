package com.hartwig.hmftools.common.variant.predicate;

import static com.hartwig.hmftools.common.variant.GermlineSampleDataTest.GERMLINE_SAMPLE_BUILDER;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.variant.GermlineSampleData;

import com.hartwig.hmftools.common.variant.GermlineSampleDataTest;
import org.junit.Test;

public class ValidGermlineSampleDataPredicateTest {

    private static final GermlineSampleData VALID = GERMLINE_SAMPLE_BUILDER.genoType("1/1").build();
    private static final GermlineSampleData INVALID = GERMLINE_SAMPLE_BUILDER.genoType("./.").build();

    @Test
    public void canCheckGermlineSampleDataCorrectly() {
        final Predicate<GermlineSampleData> predicate = new ValidGermlineSampleDataPredicate();
        assertTrue(predicate.test(VALID));
        assertFalse(predicate.test(INVALID));
    }
}
