package com.hartwig.hmftools.common.variant.predicate;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantAnnotation;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class PassFilterPredicateTest {

    private final Predicate<SomaticVariant> predicate = new PassFilterPredicate<>();

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

    private static class TestVariant implements SomaticVariant {

        @NotNull
        private final String filter;

        private TestVariant(@NotNull final String filter) {
            this.filter = filter;
        }

        @NotNull
        @Override
        public String chromosome() {
            return "1";
        }

        @Override
        public long position() {
            return 0;
        }

        @NotNull
        @Override
        public String ref() {
            return "A";
        }

        @NotNull
        @Override
        public String alt() {
            return "T";
        }

        @NotNull
        @Override
        public VariantType type() {
            return VariantType.SNP;
        }

        @NotNull
        @Override
        public String filter() {
            return filter;
        }

        @Nullable
        @Override
        public String dbsnpID() {
            return null;
        }

        @Nullable
        @Override
        public String cosmicID() {
            return null;
        }

        @NotNull
        @Override
        public List<VariantAnnotation> annotations() {
            return Lists.newArrayList();
        }

        @NotNull
        @Override
        public List<String> callers() {
            return Lists.newArrayList();
        }

        @Override
        public int totalReadCount() {
            return 0;
        }

        @Override
        public int alleleReadCount() {
            return 0;
        }
    }
}
