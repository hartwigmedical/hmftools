package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class MutationalLoadAnalyzer {

    static int analyzeVariants(@NotNull final List<SomaticVariant> variants) {
        int mutationalLoadSize = 0;
        for (final SomaticVariant variant : variants) {
            if (mutationalLoadCheck(variant)) {
                mutationalLoadSize++;
            }
        }
        return mutationalLoadSize;
    }

    private static boolean mutationalLoadCheck(@NotNull final SomaticVariant variant) {
        // KODU: Patient reporting pre-filters already, so below assert should always hold.
        assert variant.filter().equals("PASS");
        return variant.worstCodingEffect() == CodingEffect.MISSENSE;
    }
}
