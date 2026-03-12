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
    public static final double RANGE_MIN = 1;
    public static final double RANGE_MAX = 120;
    // TODO: Lookup proper threshold
    public static final double THRESHOLD = 10;

    public enum Status
    {
        HIGH,
        LOW
    }
}
