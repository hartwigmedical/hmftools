package com.hartwig.hmftools.common.variant.predicate;

import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.GermlineSampleData;
import com.hartwig.hmftools.common.variant.Variant;

import org.jetbrains.annotations.NotNull;

public final class VariantFilter {

    private VariantFilter() {
    }

    public static <T extends Variant> boolean isPass(@NotNull final T variant) {
        return new PassFilterPredicate<T>().test(variant);
    }

    @NotNull
    public static <T extends Variant> List<T> passOnly(@NotNull final List<T> variants) {
        return filter(variants, new PassFilterPredicate<>());
    }

    public static boolean validGermlineData(@NotNull final GermlineSampleData sampleData) {
        return new ValidGermlineSampleDataPredicate().test(sampleData);
    }

    @NotNull
    public static <T extends Variant> List<T> filter(@NotNull final List<T> variants,
            @NotNull final Predicate<T> predicate) {
        return variants.stream().filter(predicate).collect(Collectors.toList());
    }
}
