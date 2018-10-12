package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public final class MutationalBurdenAnalyzer {

    private static final double NUMBER_OF_MB_PER_GENOME = 2859D;

    private MutationalBurdenAnalyzer() {
    }

    static double determineTumorMutationalBurden(@NotNull List<? extends SomaticVariant> variants) {
        int tumorMutationalBurden = 0;
        for (final SomaticVariant variant : variants) {
            if (countsPassVariants(variant)) {
                tumorMutationalBurden++;
            }
        }
        return (double) tumorMutationalBurden / NUMBER_OF_MB_PER_GENOME;
    }

    private static boolean countsPassVariants(@NotNull SomaticVariant variant) {
        return !variant.isFiltered();
    }
}
