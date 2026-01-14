package com.hartwig.hmftools.datamodel.finding;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public record FindingItem<T>(
        @NotNull FindingsStatus status,
        @Nullable T finding
) {}
