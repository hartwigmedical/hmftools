package com.hartwig.hmftools.finding.datamodel;

import jakarta.validation.constraints.NotNull;

public record VisualisationFile(
        @NotNull String fileName
)
{
}
