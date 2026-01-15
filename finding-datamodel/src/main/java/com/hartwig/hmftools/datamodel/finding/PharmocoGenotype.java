package com.hartwig.hmftools.datamodel.finding;

import org.jetbrains.annotations.NotNull;

import io.soabase.recordbuilder.core.RecordBuilder;

@RecordBuilder
public record PharmocoGenotype(
        @NotNull String findingKey,
        @NotNull String gene,
        @NotNull String allele,
        int alleleCount,
        @NotNull String function,
        @Deprecated @NotNull String haplotype,
        @NotNull String linkedDrugs,
        @NotNull String urlPrescriptionInfo
) implements Finding
{
}
