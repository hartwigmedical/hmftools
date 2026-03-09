package com.hartwig.hmftools.finding.datamodel;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record PredictedTumorOrigin(
        @NotNull String findingKey,
        @NotNull CuppaMode mode,
        @NotNull String cancerType,
        double likelihood,
        @Nullable Double snvPairwiseClassifier,
        @Nullable Double genomicPositionClassifier,
        @Nullable Double featureClassifier,
        @Nullable Double altSjCohortClassifier,
        @Nullable Double expressionPairwiseClassifier
) implements Finding
{
    public enum CuppaMode
    {
        WGS,
        WGTS
    }
}
