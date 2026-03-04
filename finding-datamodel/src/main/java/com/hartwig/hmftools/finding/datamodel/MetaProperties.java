package com.hartwig.hmftools.finding.datamodel;

import java.time.LocalDate;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record MetaProperties(
        @NotNull SequencingScope sequencingScope,
        @Nullable String pipelineVersion,
        @NotNull String version,
        @NotNull RefGenomeVersion refGenomeVersion,
        @NotNull String sampleId,
        @NotNull LocalDate samplingDate)
{
}
