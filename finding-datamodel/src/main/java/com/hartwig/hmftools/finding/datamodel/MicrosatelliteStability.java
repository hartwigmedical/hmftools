package com.hartwig.hmftools.finding.datamodel;

import java.util.List;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record MicrosatelliteStability(
        @NotNull String findingKey,
        double microsatelliteIndelsPerMb,
        @NotNull MicrosatelliteStatus microsatelliteStatus,
        @NotNull List<GainDeletion> lohCopyNumbers,
        @NotNull List<String> relatedGenes
) implements Finding
{
    public enum MicrosatelliteStatus
    {
        MSI("Unstable"),
        MSS("Stable"),
        UNKNOWN("Unknown");

        @NotNull
        private final String display;

        MicrosatelliteStatus(@NotNull final String display)
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
