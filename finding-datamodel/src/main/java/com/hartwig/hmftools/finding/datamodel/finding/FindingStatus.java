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

    public enum ResultIssue
    {
        REF_REQUIRED,
        WGS_REQUIRED,
        DELETED_GENES,
        HIGH_COPY_NUMBER_NOISE,
        GENDER_MISMATCH,
        LOW_PURITY,
        TUMOR_IN_NORMAL_CONTAMINATION,
        CONTAMINATION,
        NO_TUMOR
    }

    public enum ResultStatus
    {
        NOT_AVAILABLE,
        NOT_RELIABLE,
        OK
    }
}
