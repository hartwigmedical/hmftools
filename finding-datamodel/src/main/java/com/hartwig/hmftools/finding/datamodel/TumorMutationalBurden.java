package com.hartwig.hmftools.finding.datamodel;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record TumorMutationalBurden(
        @NotNull String findingKey,
        @NotNull Status status,
        double burdenPerMb,
        int svBurden
) implements Finding
{
    public enum Status
    {
        HIGH,
        LOW
    }
}
