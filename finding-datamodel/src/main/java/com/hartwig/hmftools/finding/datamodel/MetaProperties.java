package com.hartwig.hmftools.finding.datamodel;

import java.time.LocalDate;
import java.util.SortedSet;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record MetaProperties(
        @NotNull SequencingScope sequencingScope,
        @Nullable String pipelineVersion,
        @NotNull RefGenomeVersion refGenomeVersion,
        @NotNull String sampleId,
        @NotNull LocalDate samplingDate,
        @NotNull SortedSet<String> potentialHRDGenes,
        @NotNull SortedSet<String> potentialMSIGenes)
{
}
