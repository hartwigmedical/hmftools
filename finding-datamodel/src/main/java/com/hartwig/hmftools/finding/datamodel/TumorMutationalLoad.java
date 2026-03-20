package com.hartwig.hmftools.finding.datamodel;

import com.hartwig.hmftools.finding.datamodel.finding.Finding;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record TumorMutationalLoad(
        @NotNull String findingKey,
        @NotNull Status status,
        @NotNull ThresholdValue load
) implements Finding
{

    public enum Status
    {
        HIGH,
        LOW
    }
}
