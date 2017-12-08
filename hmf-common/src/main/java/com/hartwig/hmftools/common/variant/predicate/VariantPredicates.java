package com.hartwig.hmftools.common.variant.predicate;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public final class VariantPredicates {

    private VariantPredicates() {
    }

    @SafeVarargs
    @NotNull
    public static Predicate<SomaticVariant> and(@NotNull final Predicate<SomaticVariant>... predicates) {
        return variant -> {
            for (final Predicate<SomaticVariant> predicate : predicates) {
                if (!predicate.test(variant)) {
                    return false;
                }
            }
            return true;
        };
    }

    @NotNull
    public static Predicate<SomaticVariant> not(@NotNull final Predicate<SomaticVariant> predicate) {
        return variant -> !predicate.test(variant);
    }

    @NotNull
    public static Predicate<SomaticVariant> withType(@NotNull final VariantType type) {
        return variant -> variant.type().equals(type);
    }

    @NotNull
    public static Predicate<SomaticVariant> inDBSNPAndNotInCOSMIC() {
        return and(inDBSNP(), not(inCOSMIC()));
    }

    @NotNull
    private static Predicate<SomaticVariant> inCOSMIC() {
        return SomaticVariant::isCOSMIC;
    }

    @NotNull
    private static Predicate<SomaticVariant> inDBSNP() {
        return SomaticVariant::isDBSNP;
    }
}
