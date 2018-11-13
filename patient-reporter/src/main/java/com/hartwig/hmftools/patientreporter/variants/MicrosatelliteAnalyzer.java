package com.hartwig.hmftools.patientreporter.variants;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public final class MicrosatelliteAnalyzer {

    private static final double NUMBER_OF_MB_PER_GENOME = 2859D;

    private static final int MIN_SEQUENCE_LENGTH_FOR_LONG_REPEATS = 2;
    private static final int MAX_SEQUENCE_LENGTH_FOR_LONG_REPEATS = 4;
    private static final int MIN_REPEAT_COUNT_FOR_LONG_REPEATS = 4;

    private static final int MIN_REPEAT_COUNT_FOR_SHORT_REPEATS = 5;

    private MicrosatelliteAnalyzer() {
    }

    static double determineMicrosatelliteIndelsPerMb(@NotNull List<EnrichedSomaticVariant> variants) {
        int indelCount = 0;
        for (EnrichedSomaticVariant variant : variants) {
            if (isPassIndel(variant) && repeatContextIsRelevant(variant.repeatCount(), variant.repeatSequence())) {
                indelCount++;
            }
        }
        return (double) indelCount / NUMBER_OF_MB_PER_GENOME;
    }

    private static boolean isPassIndel(@NotNull SomaticVariant variant) {
        return variant.type() == VariantType.INDEL && variant.ref().length() < 50 && variant.alt().length() < 50 && !variant.isFiltered();
    }

    @VisibleForTesting
    static boolean repeatContextIsRelevant(int repeatCount, @NotNull String sequence) {
        final int repeatSequenceLength = sequence.length();
        final boolean longRepeatRelevant =
                repeatSequenceLength >= MIN_SEQUENCE_LENGTH_FOR_LONG_REPEATS && repeatSequenceLength <= MAX_SEQUENCE_LENGTH_FOR_LONG_REPEATS
                        && repeatCount >= MIN_REPEAT_COUNT_FOR_LONG_REPEATS;
        final boolean shortRepeatRelevant = repeatSequenceLength == 1 && repeatCount >= MIN_REPEAT_COUNT_FOR_SHORT_REPEATS;
        return longRepeatRelevant | shortRepeatRelevant;
    }
}
