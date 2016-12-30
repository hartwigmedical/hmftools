package com.hartwig.hmftools.common.variant.predicate;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.variant.GermlineVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public class VCFGermlineVariantPredicate implements Predicate<GermlineVariant> {

    private static final String INVALID_VARIANT_1 = "./.";
    private static final String INVALID_VARIANT_2 = "0/0";

    @NotNull
    private final VariantType variantType;
    private final boolean isRefSample;

    public VCFGermlineVariantPredicate(@NotNull final VariantType variantType, final boolean isRefSample) {
        this.variantType = variantType;
        this.isRefSample = isRefSample;
    }

    @Override
    public boolean test(@NotNull final GermlineVariant germlineVariant) {
        boolean isVariant = false;
        if (germlineVariant.type() == variantType) {
            String dataToCheck = germlineVariant.refData();

            if (!isRefSample) {
                dataToCheck = germlineVariant.tumorData();
            }

            if (!dataToCheck.startsWith(INVALID_VARIANT_1) && !dataToCheck.startsWith(INVALID_VARIANT_2)) {
                isVariant = true;
            }
        }
        return isVariant;
    }
}
