package com.hartwig.hmftools.common.variant.predicate;

import java.util.function.Predicate;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

class PassFilterPredicate<T extends SomaticVariant> implements Predicate<T> {

    @VisibleForTesting
    static final String PASS_IDENTIFIER_1 = "PASS";
    @VisibleForTesting
    static final String PASS_IDENTIFIER_2 = ".";

    @Override
    public boolean test(@NotNull final SomaticVariant variant) {
        return variant.filter().equals(PASS_IDENTIFIER_1) || variant.filter().equals(PASS_IDENTIFIER_2);
    }
}
