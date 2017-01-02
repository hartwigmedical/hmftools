package com.hartwig.hmftools.common.variant.predicate;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public class ValidGermlineVariantPredicate implements Predicate<GermlineVariant> {

    private static final String INVALID_VARIANT_1 = "./.";
    private static final String INVALID_VARIANT_2 = "0/0";

    @NotNull
    private final VariantType variantType;
    private final boolean isRefSample;

    public ValidGermlineVariantPredicate(@NotNull final VariantType variantType, final boolean isRefSample) {
        this.variantType = variantType;
        this.isRefSample = isRefSample;
    }

    @Override
    public boolean test(@NotNull final GermlineVariant germlineVariant) {
        boolean isValid = false;
        if (germlineVariant.type() == variantType) {
            final String dataToCheck = isRefSample ? germlineVariant.refData() : germlineVariant.tumorData();

            if (!dataToCheck.startsWith(INVALID_VARIANT_1) && !dataToCheck.startsWith(INVALID_VARIANT_2)) {
                isValid = true;
            }
        }
        return isValid;
    }
}
