package com.hartwig.hmftools.finding.datamodel;

import java.util.List;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record MicrosatelliteStability(
        @NotNull String findingKey,
        @NotNull MicrosatelliteStability.Status status,
        double indelsPerMb,
        @NotNull List<GainDeletion> lohCopyNumbers,
        @NotNull List<String> relatedGenes
) implements Finding
{
    public enum Status
    {
        MSI,
        MSS
    }
}
