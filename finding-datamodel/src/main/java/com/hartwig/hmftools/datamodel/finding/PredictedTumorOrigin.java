package com.hartwig.hmftools.datamodel.finding;

import org.jetbrains.annotations.NotNull;

import io.soabase.recordbuilder.core.RecordBuilder;

@RecordBuilder
public record PredictedTumorOrigin(
        @NotNull String findingKey,
        @NotNull String cancerType,
        double likelihood
) implements Finding {}
