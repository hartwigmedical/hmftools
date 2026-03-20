package com.hartwig.hmftools.finding.datamodel;

import com.hartwig.hmftools.finding.datamodel.finding.Finding;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record TumorMutationalBurden(
        @NotNull String findingKey,
        @NotNull Status status,
        @NotNull ThresholdValue burdenPerMb,
        int svBurden
) implements Finding
{
    public enum Status
    {
        HIGH,
        LOW
    }
}
