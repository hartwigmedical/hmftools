package com.hartwig.hmftools.common.variant.predicate;

import java.util.function.Predicate;

import com.hartwig.hmftools.common.variant.GermlineSampleData;

import org.jetbrains.annotations.Nullable;

class ValidGermlineSampleDataPredicate implements Predicate<GermlineSampleData> {

    private static final String INVALID_GENOTYPE_1 = "./.";
    private static final String INVALID_GENOTYPE_2 = "0/0";

    @Override
    public boolean test(@Nullable final GermlineSampleData data) {
        return (data != null && !data.genoType().equals(INVALID_GENOTYPE_1) && !data.genoType().equals(
                INVALID_GENOTYPE_2));
    }
}
