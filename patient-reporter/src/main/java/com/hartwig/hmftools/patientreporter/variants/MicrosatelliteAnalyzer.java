package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public final class MicrosatelliteAnalyzer {

    private MicrosatelliteAnalyzer() {
    }

    static double determineMicrosatelliteIndels(@NotNull final List<EnrichedSomaticVariant> variants) {
        int indelCount = 0;
        for (EnrichedSomaticVariant variant : variants) {
            if (isPassIndel(variant) && repeatContextIsRelevant(variant.repeatCount(), variant.repeatSequence())) {
                indelCount++;
            }
        }
        return (double) indelCount / 3095D;
    }

    private static boolean isPassIndel(@NotNull final SomaticVariant variant) {
        return variant.type() == VariantType.INDEL && variant.ref().length() < 50 && variant.alt().length() < 50 && !variant.isFiltered();
    }

    @VisibleForTesting
    static boolean repeatContextIsRelevant(int repeatCount, @NotNull String sequence) {
        final int repeatSequenceLength = sequence.length();
        final boolean longRepeatRelevant = repeatSequenceLength >= 2 && repeatSequenceLength <= 4 && repeatCount >= 4;
        final boolean shortRepeatRelevant = repeatSequenceLength == 1 && repeatCount >= 5;
        return longRepeatRelevant | shortRepeatRelevant;
    }
}
