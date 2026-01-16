package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;

import org.jspecify.annotations.Nullable;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record MetaProperties(
        @NotNull ExperimentType experimentType,
        @Nullable String pipelineVersion,
        @NotNull String version,
        @NotNull OrangeRefGenomeVersion refGenomeVersion
)
{
}
