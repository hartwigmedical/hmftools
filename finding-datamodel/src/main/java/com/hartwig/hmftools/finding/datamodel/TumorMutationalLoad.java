package com.hartwig.hmftools.finding.datamodel;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record TumorMutationalLoad(
        @NotNull String findingKey,
        @NotNull Status status,
        int load
) implements Finding
{
    public static final double RANGE_MIN = 1;
    public static final double RANGE_MAX = 1000;
    // TODO: Lookup proper threshold
    public static final double THRESHOLD = 10;

    public enum Status
    {
        HIGH,
        LOW
    }
}
