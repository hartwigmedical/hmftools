package com.hartwig.hmftools.common.cuppa;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class CuppaTestFactory {

    private CuppaTestFactory() { }

    @NotNull
    public static CuppaData createMinimalCuppaData()
    {
        // Some downstream algo's expect at least one prediction, so that is considered "minimal"
        List<CuppaPrediction> predictions = Lists.newArrayList();
        predictions.add(ImmutableCuppaPrediction.builder().cancerType("cancer").likelihood(1D).build());

        return ImmutableCuppaData.builder()
                .predictions(predictions)
                .simpleDups32To200B(0)
                .maxComplexSize(0)
                .telomericSGLs(0)
                .LINECount(0)
                .build();
    }
}
