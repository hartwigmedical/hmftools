package com.hartwig.hmftools.patientreporter.variants;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class MutationalLoadAnalyzer {

    static int analyzeVariants(@NotNull final List<SomaticVariant> variants) {
        final List<Boolean> mutationalLoad = new ArrayList<>();
        int mutationalLoadSize = 0;
        for (final SomaticVariant variant : variants) {
            if (mutationalLoadCheck(variant)) {
                mutationalLoad.add(mutationalLoadCheck(variant));
            }
            mutationalLoadSize++;
        }

        return mutationalLoad.size();
    }

    static boolean mutationalLoadCheck(@NotNull final SomaticVariant variant) {
        return variant.worstCodingEffect().contains("MISSENSE") && variant.filter().equals("PASS");
    }
}
