package com.hartwig.hmftools.datamodel.finding;

import org.jspecify.annotations.Nullable;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record FindingItem<T>(
        @NotNull FindingsStatus status,
        @Nullable T finding
)
{
}
