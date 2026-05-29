package com.hartwig.hmftools.finding.datamodel;

import java.util.Set;

import com.hartwig.hmftools.finding.datamodel.finding.Finding;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record RnaStatistics(
        @NotNull String findingKey,
        @NotNull Set<QcStatus> qcStatus,
        long totalFragments,
        long duplicateFragments,
        double splicedFragmentPercent,
        double unsplicedFragmentPercent,
        double altFragmentPercent,
        double chimericFragmentPercent
) implements Finding
{
    public enum QcStatus
    {
        PASS,
        FAIL_LOW_COVERAGE,
        WARN_DUPLICATE_RATE,
        WARN_SPLICED_GENE_COVERAGE
    }
}
