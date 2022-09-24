package com.hartwig.hmftools.common.metrics;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class WGSMetrics
{
    public abstract double meanCoverage();

    public abstract double sdCoverage();

    public abstract int medianCoverage();

    public abstract int madCoverage();

    @Nullable
    public abstract Double pctExcAdapter();

    public abstract double pctExcMapQ();

    public abstract double pctExcDupe();

    public abstract double pctExcUnpaired();

    public abstract double pctExcBaseQ();

    public abstract double pctExcOverlap();

    public abstract double pctExcCapped();

    public abstract double pctExcTotal();

    @Nullable
    public abstract Double coverage1xPercentage();

    public abstract double coverage10xPercentage();

    public abstract double coverage20xPercentage();

    public abstract double coverage30xPercentage();

    public abstract double coverage60xPercentage();
}
