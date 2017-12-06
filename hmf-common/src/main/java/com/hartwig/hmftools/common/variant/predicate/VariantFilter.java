package com.hartwig.hmftools.common.variant.predicate;

import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.Variant;

import org.jetbrains.annotations.NotNull;

public final class VariantFilter {

    private VariantFilter() {
    }

    @NotNull
    public static <T extends Variant> List<T> passOnly(@NotNull final List<T> variants) {
        return filter(variants, new PassFilterPredicate<>());
    }

    @NotNull
    public static <T extends Variant> List<T> filter(@NotNull final List<T> variants,
            @NotNull final Predicate<T> predicate) {
        return variants.stream().filter(predicate).collect(Collectors.toList());
    }
}
