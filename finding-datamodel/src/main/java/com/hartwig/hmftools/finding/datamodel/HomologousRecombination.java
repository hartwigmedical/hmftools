package com.hartwig.hmftools.finding.datamodel;

import java.util.List;

import com.hartwig.hmftools.finding.datamodel.finding.Finding;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record HomologousRecombination(
        @NotNull String findingKey,
        @NotNull HomologousRecombination.Status status,
        @NotNull ThresholdValue hrdValue,
        double brca1Value,
        double brca2Value,
        @NotNull String hrdType,
        @NotNull List<GainDeletion> lohCopyNumbers,
        @NotNull List<String> drivingGenes
) implements Finding
{

    public enum Status
    {
        UNDETERMINED,
        HR_PROFICIENT,
        HR_DEFICIENT
    }
}
