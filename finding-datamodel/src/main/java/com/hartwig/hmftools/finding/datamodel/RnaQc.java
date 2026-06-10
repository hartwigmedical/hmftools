package com.hartwig.hmftools.finding.datamodel;

import java.util.SortedSet;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record RnaQc(
        @NotNull SortedSet<QcStatus> warnings,
        @NotNull SortedSet<QcStatus> errors,
        long totalFragments,
        long duplicateFragments,
        double splicedFragmentPercent,
        double unsplicedFragmentPercent,
        double altFragmentPercent,
        double chimericFragmentPercent
)
{
    public enum QcStatus
    {
        WARN_DUPLICATE_RATE,
        WARN_SPLICED_GENE_COVERAGE,
        FAIL_LOW_COVERAGE
    }
}
