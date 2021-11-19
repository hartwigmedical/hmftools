package com.hartwig.hmftools.common.cuppa;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class CuppaTestFactory {

    private CuppaTestFactory() {
    }

    @NotNull
    public static CuppaData createMinimalCuppaData() {
        return ImmutableCuppaData.builder()
                .predictedCancerType(Strings.EMPTY)
                .bestPredictionLikelihood(0D)
                .simpleDups32To200B(0)
                .maxComplexSize(0)
                .telomericSGLs(0)
                .LINECount(0)
                .build();
    }
}
