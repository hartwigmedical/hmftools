package com.hartwig.hmftools.finding.datamodel;

import com.hartwig.hmftools.datamodel.cuppa.CuppaMode;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record PredictedTumorOrigin(
        @NotNull String findingKey,
        @NotNull CuppaMode mode,
        @NotNull String cancerType,
        double likelihood,
        Double snvPairwiseClassifier,
        Double genomicPositionClassifier,
        Double featureClassifier,
        Double altSjCohortClassifier,
        Double expressionPairwiseClassifier
) implements Finding
{
}
