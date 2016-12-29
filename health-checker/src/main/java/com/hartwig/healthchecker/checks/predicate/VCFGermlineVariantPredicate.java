package com.hartwig.healthchecker.checks.predicate;

import java.util.function.Predicate;

import com.hartwig.healthchecker.checks.model.VCFGermlineData;
import com.hartwig.healthchecker.checks.model.VCFType;

import org.jetbrains.annotations.NotNull;

public class VCFGermlineVariantPredicate implements Predicate<VCFGermlineData> {

    private static final String INVALID_VARIANT_1 = "./.";
    private static final String INVALID_VARIANT_2 = "0/0";

    @NotNull
    private final VCFType vcfType;
    private final boolean isRefSample;

    public VCFGermlineVariantPredicate(@NotNull final VCFType vcfType, final boolean isRefSample) {
        this.vcfType = vcfType;
        this.isRefSample = isRefSample;
    }

    @Override
    public boolean test(@NotNull final VCFGermlineData vcfGermlineData) {
        boolean isVariant = false;
        if (vcfGermlineData.getType() == vcfType) {
            String dataToCheck = vcfGermlineData.getRefData();

            if (!isRefSample) {
                dataToCheck = vcfGermlineData.getTumorData();
            }

            if (!dataToCheck.startsWith(INVALID_VARIANT_1) && !dataToCheck.startsWith(INVALID_VARIANT_2)) {
                isVariant = true;
            }
        }
        return isVariant;
    }
}
