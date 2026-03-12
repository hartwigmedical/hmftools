package com.hartwig.hmftools.finding.datamodel;

import java.util.List;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record HomologousRecombination(
        @NotNull String findingKey,
        double brca1Value,
        double brca2Value,
        double hrdValue,
        @NotNull HomologousRecombination.Status status,
        @NotNull String hrdType,
        @NotNull List<GainDeletion> lohCopyNumbers,
        @NotNull List<String> relatedGenes
) implements Finding
{
    public enum Status
    {
        HR_PROFICIENT,
        HR_DEFICIENT
    }
}
