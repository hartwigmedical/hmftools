package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

import org.jetbrains.annotations.NotNull;

import io.soabase.recordbuilder.core.RecordBuilder;

@RecordBuilder
public record FindingList<T extends Finding>(
        @NotNull FindingsStatus status,
        @NotNull List<T> all)
{
}
