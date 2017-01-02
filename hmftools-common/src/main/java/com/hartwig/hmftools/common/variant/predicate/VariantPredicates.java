package com.hartwig.hmftools.common.variant.predicate;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public final class VariantPredicates {

    private VariantPredicates() {
    }

    @NotNull
    public static Predicate<SomaticVariant> withType(@NotNull final VariantType type) {
        return variant -> variant.type().equals(type);
    }

    @NotNull
    public static Predicate<SomaticVariant> inDBSNPAndNotInCOSMIC() {
        return variant -> variant.isDBSNP() && !variant.isCOSMIC();
    }

    @NotNull
    public static Predicate<SomaticVariant> withCaller(@NotNull final String caller) {
        return variant -> variant.callers().contains(caller);
    }

    @NotNull
    public static Predicate<SomaticVariant> withMinCallers(final int count) {
        return variant -> variant.callerCount() >= count;
    }

    @NotNull
    public static Predicate<SomaticVariant> withExactCallerCount(final int count) {
        return variant -> variant.callerCount() == count;
    }
}
