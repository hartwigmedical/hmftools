package com.hartwig.hmftools.finding.datamodel;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record TumorMutationStatus(
        @NotNull String findingKey,
        double tumorMutationalBurdenPerMb,
        @NotNull Status tumorMutationalBurdenStatus,
        int tumorMutationalLoad,
        @NotNull Status tumorMutationalLoadStatus,
        int svTumorMutationalBurden
) implements Finding
{
    public enum Status
    {
        HIGH("High"),
        LOW("Low"),
        UNKNOWN("Unknown");

        @NotNull
        private final String display;

        Status(@NotNull final String display)
        {
            this.display = display;
        }

        @NotNull
        public String display()
        {
            return display;
        }
    }
}
