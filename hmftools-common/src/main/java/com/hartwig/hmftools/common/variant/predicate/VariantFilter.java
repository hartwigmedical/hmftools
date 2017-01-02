package com.hartwig.hmftools.common.variant.predicate;

import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public final class VariantFilter {

    private VariantFilter() {
    }

    @NotNull
    public static <T extends Variant> List<T> passOnly(@NotNull List<T> variants) {
        return filter(variants, new PassFilterPredicate<>());
    }

    @NotNull
    public static List<GermlineVariant> validGermline(@NotNull List<GermlineVariant> variants,
            @NotNull VariantType type, boolean isRefSample) {
        return filter(variants, new ValidGermlineVariantPredicate(type, isRefSample));
    }

    @NotNull
    public static <T extends Variant> List<T> filter(@NotNull List<T> variants, @NotNull Predicate<T> predicate) {
        return variants.stream().filter(predicate).collect(Collectors.toList());
    }
}
