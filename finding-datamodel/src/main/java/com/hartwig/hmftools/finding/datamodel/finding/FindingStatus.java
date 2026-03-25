package com.hartwig.hmftools.finding.datamodel.finding;

import java.util.SortedSet;

import com.hartwig.hmftools.finding.datamodel.RecordBuilder;

import jakarta.validation.constraints.NotNull;

@SuppressWarnings("unused")
@RecordBuilder
public record FindingStatus(@NotNull FindingStatus.Status status, @NotNull SortedSet<Issue> errors, @NotNull SortedSet<Issue> warnings)
{
    public boolean isOK()
    {
        return status == Status.OK;
    }

    public boolean isNotAvailable()
    {
        return status == Status.NOT_AVAILABLE;
    }

    public boolean isNotReliable()
    {
        return status == Status.NOT_RELIABLE;
    }

    public enum Issue
    {
        REF_REQUIRED,
        WGS_REQUIRED,
        DELETED_GENES,
        HIGH_COPY_NUMBER_NOISE,
        GENDER_MISMATCH,
        LOW_PURITY,
        TUMOR_IN_NORMAL_CONTAMINATION,
        CONTAMINATION,
        NO_TUMOR,
        // No predicted tumor origins that meet likelihood threshold
        NO_REPORTABLE_VALUE
    }

    public enum Status
    {
        NOT_AVAILABLE,
        NOT_RELIABLE,
        OK
    }
}
