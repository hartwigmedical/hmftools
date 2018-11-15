package com.hartwig.hmftools.strelka.mnv.scores;

import com.hartwig.hmftools.common.sam.SamRecords;

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
        return ImmutableVariantScore.of(type, SamRecords.avgQuality(baseQualities));
    }
}
