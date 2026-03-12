package com.hartwig.hmftools.finding.datamodel;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record TumorMutationalLoad(
        @NotNull String findingKey,
        @NotNull Status status,
        int load
) implements Finding
{
    public enum Status
    {
        HIGH,
        LOW
    }
}
