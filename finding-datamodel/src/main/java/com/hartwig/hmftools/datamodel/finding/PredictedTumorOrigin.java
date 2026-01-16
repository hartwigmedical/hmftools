package com.hartwig.hmftools.datamodel.finding;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record PredictedTumorOrigin(
        @NotNull String findingKey,
        @NotNull String cancerType,
        double likelihood
) implements Finding
{
}
