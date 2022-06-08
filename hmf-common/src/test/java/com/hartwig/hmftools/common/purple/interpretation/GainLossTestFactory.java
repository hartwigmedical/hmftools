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
        return ImmutableGainLoss.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(gene)
                .transcript(Strings.EMPTY)
                .isCanonical(true)
                .interpretation(interpretation)
                .minCopies(1)
                .maxCopies(1)
                .build();
    }
}
