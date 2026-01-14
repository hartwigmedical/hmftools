package com.hartwig.hmftools.datamodel.finding;

import org.jetbrains.annotations.NotNull;

public record PredictedTumorOrigin(
        @NotNull String findingKey,
        @NotNull String cancerType,
        double likelihood
) implements Finding {}
