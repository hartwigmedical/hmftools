package com.hartwig.hmftools.common.metrics;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class WGSMetricWithQC
{

    @NotNull
    public abstract WGSMetrics refMetrics();

    @NotNull
    public abstract WGSMetrics tumorMetrics();

    public abstract boolean sufficientCoverage();
}
