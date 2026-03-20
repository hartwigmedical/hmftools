package com.hartwig.hmftools.finding.datamodel;

import java.util.List;

import com.hartwig.hmftools.finding.datamodel.finding.Finding;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record PredictedTumorOrigin(
        @NotNull String findingKey,
        @NotNull CuppaMode mode,
        @NotNull List<Prediction> predictions,
        @Nullable VisualisationFile visualisationFile
) implements Finding
{
    public enum CuppaMode
    {
        WGS,
        WGTS
    }

    @RecordBuilder
    public record Prediction(
            @NotNull String findingKey,
            @NotNull String cancerType,
            double likelihood,
            @Nullable Double snvPairwiseClassifier,
            @Nullable Double genomicPositionClassifier,
            @Nullable Double featureClassifier,
            @Nullable Double altSjCohortClassifier,
            @Nullable Double expressionPairwiseClassifier
    ) implements Finding
    {
    }
}
