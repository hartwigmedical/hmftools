package com.hartwig.hmftools.finding.datamodel;

import java.util.List;

import com.hartwig.hmftools.finding.datamodel.finding.Finding;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record MicrosatelliteStability(
        @NotNull String findingKey,
        @NotNull MicrosatelliteStability.Status status,
        @NotNull ThresholdValue indelsPerMb,
        @NotNull List<GainDeletion> lohCopyNumbers,
        @NotNull List<String> drivingGenes
) implements Finding
{

    public enum Status
    {
        MSI,
        MSS
    }
}
