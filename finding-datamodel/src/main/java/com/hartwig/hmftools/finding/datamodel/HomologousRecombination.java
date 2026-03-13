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
        @NotNull List<String> drivingGenes
) implements Finding
{
    public static final double RANGE_MIN = 0;
    public static final double RANGE_MAX = 1;
    // TODO: Lookup proper threshold
    public static final double THRESHOLD = 0.5;

    public enum Status
    {
        UNDETERMINED,
        HR_PROFICIENT,
        HR_DEFICIENT
    }
}
