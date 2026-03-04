package com.hartwig.hmftools.finding.datamodel;

import java.util.List;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record HomologousRecombination(
        @NotNull String findingKey,
        double brca1Value,
        double brca2Value,
        double hrdValue,
        @NotNull HrStatus hrStatus,
        @NotNull String hrdType,
        @NotNull List<GainDeletion> lohCopyNumbers,
        @NotNull List<String> genes
) implements Finding
{
    public enum HrStatus
    {
        CANNOT_BE_DETERMINED("Cannot be determined"),
        HR_PROFICIENT("Proficient"),
        HR_DEFICIENT("Deficient"),
        UNKNOWN("Unknown");

        @NotNull
        private final String display;

        HrStatus(@NotNull final String display)
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
