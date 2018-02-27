package com.hartwig.hmftools.common.variant.predicate;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.ImmutableSomaticVariantImpl;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantImpl;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.snpeff.VariantAnnotation;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class PassFilterPredicateTest {

    private final Predicate<SomaticVariant> predicate = new PassFilterPredicate<>();

    @Test
    public void filtersOnFirstPassIdentifier() {
        assertTrue(predicate.test(create(PassFilterPredicate.PASS_IDENTIFIER_1)));
    }

    @Test
    public void filtersOnSecondPassIdentifier() {
        assertTrue(predicate.test(create(PassFilterPredicate.PASS_IDENTIFIER_2)));
    }

    @Test
    public void filterOnNonPassing() {
        assertFalse(predicate.test(create("NotPassing!")));
    }

    private SomaticVariant create(@NotNull final String filter) {
        return new SomaticVariantImpl.Builder().filter(filter).build();
    }

}
