package com.hartwig.hmftools.finding.datamodel.finding;

import java.util.SortedSet;

import com.hartwig.hmftools.finding.datamodel.RecordBuilder;

import jakarta.validation.constraints.NotNull;

@SuppressWarnings("unused")
@RecordBuilder
public record FindingStatus(@NotNull FindingStatus.Status status, @NotNull SortedSet<Issue> errors, @NotNull SortedSet<Issue> warnings)
{
    public enum Status
    {
        NOT_AVAILABLE,
        NOT_RELIABLE,
        OK
    }

    // TODO: Should use Qc.Status for this?
    public enum Issue
    {
        // TODO: inconsistent naming ref/normal (vs Qc.Status)
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
        NO_REPORTABLE_VALUE,
        TUMOR_SAMPLE_QUALITY_CONTROL,
        REF_SAMPLE_QUALITY_CONTROL
    }

    public boolean isOK()
    {
        return status == Status.OK;
    }

    public boolean isAvailable()
    {
        return status == Status.OK || status == Status.NOT_RELIABLE;
    }

    public boolean isReliable()
    {
        return status != Status.NOT_RELIABLE;
    }
}
