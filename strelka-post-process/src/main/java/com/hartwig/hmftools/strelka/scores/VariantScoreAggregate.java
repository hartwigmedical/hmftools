package com.hartwig.hmftools.strelka.scores;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

@Value.Immutable
@Value.Style(allParameters = true)
public abstract class VariantScoreAggregate {

    public abstract long sum();

    public abstract int count();

    @Value.Lazy
    public int average() {
        if (count() == 0) {
            return 0;
        } else {
            return (int) sum() / count();
        }
    }

    public VariantScoreAggregate addScore(@NotNull final VariantScore variantScore) {
        if (variantScore.type() == ReadType.ALT || variantScore.type() == ReadType.REF) {
            return ImmutableVariantScoreAggregate.of(sum() + variantScore.score(), count() + 1);
        } else {
            return this;
        }
    }

    public static VariantScoreAggregate newScore() {
        return ImmutableVariantScoreAggregate.of(0, 0);
    }
}
