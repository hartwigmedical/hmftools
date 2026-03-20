package com.hartwig.hmftools.finding.datamodel.finding;

import com.hartwig.hmftools.finding.datamodel.RecordBuilder;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record FindingItem<T>(
        @NotNull FindingsStatus status,
        @Nullable T finding
)
{
    @Nullable
    public T findingIfOK()
    {
        return status.isOK() ? finding : null;
    }
}
