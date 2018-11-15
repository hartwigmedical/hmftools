package com.hartwig.hmftools.strelka.mnv.scores;

import static com.hartwig.hmftools.common.sam.SAMRecords.getBaseQuality;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class VariantScore {
    public abstract ReadType type();

    public abstract int score();

    public static VariantScore of(@NotNull final ReadType type, final char baseScore) {
        return of(type, String.valueOf(baseScore));
    }

    public static VariantScore of(@NotNull final ReadType type, @NotNull final String baseQualities) {
        return ImmutableVariantScore.of(type, avgQuality(baseQualities));
    }

    private static int avgQuality(@NotNull final String baseQualities) {
        return (int) Math.floor(totalQuality(baseQualities) / baseQualities.length());
    }

    private static int totalQuality(@NotNull final String baseQualities) {
        int score = 0;
        for (int index = 0; index < baseQualities.length(); index++) {
            score += getBaseQuality(baseQualities.charAt(index));
        }
        return score;
    }
}
