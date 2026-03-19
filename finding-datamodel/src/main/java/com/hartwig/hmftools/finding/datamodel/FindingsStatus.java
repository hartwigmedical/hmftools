package com.hartwig.hmftools.finding.datamodel;

import java.util.SortedSet;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record FindingsStatus(@NotNull ResultStatus status, @NotNull SortedSet<ResultIssue> errors, @NotNull SortedSet<ResultIssue> warnings)
{
    public boolean isOK()
    {
        return status == ResultStatus.OK;
    }
}
