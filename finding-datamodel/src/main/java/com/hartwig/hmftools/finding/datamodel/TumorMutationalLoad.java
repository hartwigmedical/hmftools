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
        HIGH("High"),
        LOW("Low");

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
