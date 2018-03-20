package com.hartwig.hmftools.common.metrics;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class WGSMetrics {

    public abstract double refMeanCoverage();

    public abstract double ref10xCoveragePercentage();

    public abstract double ref20xCoveragePercentage();

    @Nullable
    public abstract Double tumorMeanCoverage();

    @Nullable
    public abstract Double tumor30xCoveragePercentage();

    @Nullable
    public abstract Double tumor60xCoveragePercentage();
}
