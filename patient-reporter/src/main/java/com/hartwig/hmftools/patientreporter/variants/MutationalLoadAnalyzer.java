package com.hartwig.hmftools.patientreporter.variants;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public class MutationalLoadAnalyzer {

    static int analyzeVariants(@NotNull final List<SomaticVariant> variants) {
        final List<Boolean> mutationalLoadList = new ArrayList<>();
        int mutationalLoadSize = 0;
        for (final SomaticVariant variant : variants) {
            final boolean mutationLoadValue = mutationalLoad(variant);
            mutationalLoadList.add(mutationLoadValue);
            mutationalLoadSize++;
        }

        return mutationalLoadList.size();
    }

    static boolean mutationalLoad(@NotNull final SomaticVariant variant) {
        return variant.worstCodingEffect().contains("missense") && variant.filter().equals("PASS");
    }
}
