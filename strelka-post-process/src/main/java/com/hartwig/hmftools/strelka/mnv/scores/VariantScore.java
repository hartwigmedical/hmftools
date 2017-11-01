package com.hartwig.hmftools.strelka.mnv.scores;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class VariantScore {
    private static final int PHRED_OFFSET = 33;

    public abstract ReadType type();

    public abstract int score();

    public static VariantScore of(@NotNull final ReadType type, final char baseScore) {
        return ImmutableVariantScore.of(type, baseScore - PHRED_OFFSET);
    }

    public static VariantScore of(@NotNull final ReadType type, @NotNull final String baseQualities) {
        double score = 0;
        for (int index = 0; index < baseQualities.length(); index++) {
            score += baseQualities.charAt(index) - PHRED_OFFSET;
        }
        return ImmutableVariantScore.of(type, (int) Math.floor(score / baseQualities.length()));
    }
}
