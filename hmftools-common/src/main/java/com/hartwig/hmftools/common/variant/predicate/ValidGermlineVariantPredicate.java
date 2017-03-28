package com.hartwig.hmftools.common.variant.predicate;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.variant.GermlineSampleData;
import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

class ValidGermlineVariantPredicate implements Predicate<GermlineVariant> {

    private static final String INVALID_GENOTYPE_1 = "./.";
    private static final String INVALID_GENOTYPE_2 = "0/0";

    @NotNull
    private final VariantType variantType;
    private final boolean isRefSample;

    ValidGermlineVariantPredicate(@NotNull final VariantType variantType, final boolean isRefSample) {
        this.variantType = variantType;
        this.isRefSample = isRefSample;
    }

    @Override
    public boolean test(@NotNull final GermlineVariant variant) {
        boolean isValid = false;
        if (variant.type() == variantType) {
            final GermlineSampleData dataToCheck = isRefSample ? variant.refData() : variant.tumorData();

            if (dataToCheck != null && !dataToCheck.genoType().equals(INVALID_GENOTYPE_1)
                    && !dataToCheck.genoType().equals(INVALID_GENOTYPE_2)) {
                isValid = true;
            }
        }
        return isValid;
    }
}
