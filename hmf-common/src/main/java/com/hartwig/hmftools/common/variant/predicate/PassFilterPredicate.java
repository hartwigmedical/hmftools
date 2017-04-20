package com.hartwig.hmftools.common.variant.predicate;

import java.util.function.Predicate;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.variant.Variant;

import org.jetbrains.annotations.NotNull;

class PassFilterPredicate<T extends Variant> implements Predicate<T> {

    @VisibleForTesting
    static final String PASS_IDENTIFIER_1 = "PASS";
    @VisibleForTesting
    static final String PASS_IDENTIFIER_2 = ".";

    @Override
    public boolean test(@NotNull final Variant variant) {
        return variant.filter().equals(PASS_IDENTIFIER_1) || variant.filter().equals(PASS_IDENTIFIER_2);
    }
}
