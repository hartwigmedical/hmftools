package com.hartwig.hmftools.finding.datamodel;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record PharmocoGenotype(
        @NotNull String findingKey,
        @NotNull String gene,
        @NotNull String allele,
        int alleleCount,
        @NotNull String function,
        @NotNull String haplotype,
        @NotNull String linkedDrugs,
        @NotNull String urlPrescriptionInfo
) implements Finding
{
}
