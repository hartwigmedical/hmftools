package com.hartwig.hmftools.common.variant.predicate;

import java.util.List;
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
        return variants.stream().filter(new PassFilterPredicate()).collect(Collectors.toList());
    }

    @NotNull
    public static List<GermlineVariant> validGermline(@NotNull List<GermlineVariant> variants,
            @NotNull VariantType type, boolean isRefSample) {
        return variants.stream().filter(new ValidGermlineVariantPredicate(type, isRefSample)).collect(
                Collectors.toList());
    }
}
