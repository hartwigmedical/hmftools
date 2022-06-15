package com.hartwig.hmftools.common.purple.interpretation;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class GainLossTestFactory {

    private GainLossTestFactory() {
    }

    @NotNull
    public static GainLoss createTestGainLoss() {
        return createGainLoss(Strings.EMPTY, CopyNumberInterpretation.FULL_GAIN);
    }

    @NotNull
    public static GainLoss createGainLoss(@NotNull String gene, @NotNull CopyNumberInterpretation interpretation) {
        return builder().gene(gene).interpretation(interpretation).build();
    }

    @NotNull
    public static ImmutableGainLoss.Builder builder() {
        return ImmutableGainLoss.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .isCanonical(true)
                .interpretation(CopyNumberInterpretation.FULL_GAIN)
                .minCopies(0)
                .maxCopies(0);
    }
}
