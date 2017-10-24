package com.hartwig.hmftools.strelka.scores;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

@Value.Immutable
@Value.Style(allParameters = true)
public abstract class VariantScore {

    public abstract long sum();

    public abstract int count();

    @Value.Lazy
    public int average() {
        return (int) sum() / count();
    }

    public VariantScore addScore(@NotNull final ReadScore readScore) {
        if (readScore.type() == ReadType.ALT || readScore.type() == ReadType.REF) {
            return ImmutableVariantScore.of(sum() + readScore.score(), count() + 1);
        } else {
            return this;
        }
    }

    public static VariantScore newScore() {
        return ImmutableVariantScore.of(0, 0);
    }
}
