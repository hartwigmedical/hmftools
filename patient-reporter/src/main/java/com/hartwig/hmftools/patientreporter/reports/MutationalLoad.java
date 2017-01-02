package com.hartwig.hmftools.patientreporter.reports;

import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.predicate.VariantFilter;

import org.jetbrains.annotations.NotNull;

public final class MutationalLoad {

    private MutationalLoad() {
    }

    public int calcMutationalLoad(@NotNull List<SomaticVariant> variants) {
        variants = VariantFilter.passOnly(variants);
        return 1;
    }
}
