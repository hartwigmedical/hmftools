package com.hartwig.hmftools.common.variant.predicate;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.variant.Variant;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PassFilterPredicateTest {

    private final Predicate<Variant> predicate = new PassFilterPredicate();

    @Test
    public void filtersOnFirstPassIdentifier() {
        assertTrue(predicate.test(new TestVariant(PassFilterPredicate.PASS_IDENTIFIER_1)));
    }

    @Test
    public void filtersOnSecondPassIdentifier() {
        assertTrue(predicate.test(new TestVariant(PassFilterPredicate.PASS_IDENTIFIER_2)));
    }

    @Test
    public void filterOnNonPassing() {
        assertFalse(predicate.test(new TestVariant("NotPassing!")));
    }

    private static class TestVariant implements Variant {

        @NotNull
        private final String filter;

        private TestVariant(@NotNull final String filter) {
            this.filter = filter;
        }

        @NotNull
        @Override
        public String filter() {
            return filter;
        }
    }
}
