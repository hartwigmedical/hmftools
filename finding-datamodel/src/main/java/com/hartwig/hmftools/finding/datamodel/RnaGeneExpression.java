package com.hartwig.hmftools.finding.datamodel;

import com.hartwig.hmftools.finding.datamodel.finding.Finding;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record RnaGeneExpression(
        @NotNull String findingKey,
        @NotNull String gene,
        double tpm,
        double medianTpmCohort,
        double percentileCohort,
        @Nullable Double medianTpmCancer,
        @Nullable Double percentileCancer
) implements Finding
{
}
