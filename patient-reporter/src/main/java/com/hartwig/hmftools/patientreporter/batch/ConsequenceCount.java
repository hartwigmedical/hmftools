package com.hartwig.hmftools.patientreporter.batch;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantConsequence;

import org.jetbrains.annotations.NotNull;

final class ConsequenceCount {

    private ConsequenceCount() {
    }

    @NotNull
    static Map<VariantConsequence, Integer> count(@NotNull final List<SomaticVariant> variants) {
        final Map<VariantConsequence, Integer> counts = Maps.newHashMap();
        for (final VariantConsequence consequence : VariantConsequence.values()) {
            counts.put(consequence, 0);
        }

        for (final SomaticVariant variant : variants) {
            for (final VariantConsequence consequence : VariantConsequence.values()) {
                if (variant.hasConsequence(consequence)) {
                    counts.put(consequence, counts.get(consequence) + 1);
                }
            }
        }
        return counts;
    }
}
