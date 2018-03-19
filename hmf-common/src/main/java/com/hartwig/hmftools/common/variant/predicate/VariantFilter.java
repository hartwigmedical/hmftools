package com.hartwig.hmftools.common.variant.predicate;

import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public final class VariantFilter {

    private VariantFilter() {
    }

    @NotNull
    public static <T extends SomaticVariant> List<T> filter(@NotNull final List<T> variants, @NotNull final Predicate<T> predicate) {
        return variants.stream().filter(predicate).collect(Collectors.toList());
    }
}
