package com.hartwig.hmftools.finding.datamodel;

import java.util.List;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record MicrosatelliteStability(
        @NotNull String findingKey,
        @NotNull MicrosatelliteStability.Status status,
        double indelsPerMb,
        @NotNull List<GainDeletion> lohCopyNumbers,
        @NotNull List<String> drivingGenes
) implements Finding
{
    public static final double RANGE_MIN = 1;
    public static final double RANGE_MAX = 100;
    // TODO: Lookup proper threshold
    public static final double THRESHOLD = 4.0;

    public enum Status
    {
        MSI,
        MSS
    }
}
