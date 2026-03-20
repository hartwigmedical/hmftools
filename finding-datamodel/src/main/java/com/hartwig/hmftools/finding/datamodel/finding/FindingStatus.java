package com.hartwig.hmftools.finding.datamodel.finding;

import java.util.SortedSet;

import com.hartwig.hmftools.finding.datamodel.RecordBuilder;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record FindingStatus(@NotNull ResultStatus status, @NotNull SortedSet<ResultIssue> errors, @NotNull SortedSet<ResultIssue> warnings)
{
    public boolean isOK()
    {
        return status == ResultStatus.OK;
    }
}
