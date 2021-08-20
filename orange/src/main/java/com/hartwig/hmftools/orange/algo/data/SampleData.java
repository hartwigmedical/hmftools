package com.hartwig.hmftools.orange.algo.data;

import com.hartwig.hmftools.common.flagstat.Flagstat;
import com.hartwig.hmftools.common.metrics.WGSMetrics;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class SampleData {

    @NotNull
    public abstract WGSMetrics metrics();

    @NotNull
    public abstract Flagstat flagstat();
}
