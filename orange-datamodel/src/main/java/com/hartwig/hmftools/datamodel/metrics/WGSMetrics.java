package com.hartwig.hmftools.datamodel.metrics;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface WGSMetrics
{
    double meanCoverage();

    double sdCoverage();

    int medianCoverage();

    int madCoverage();

    @Nullable
    Double pctExcAdapter();

    double pctExcMapQ();

    double pctExcDupe();

    double pctExcUnpaired(); // this is actually unmapped not unpaired percent

    double pctExcBaseQ();

    double pctExcOverlap();

    double pctExcCapped();

    double pctExcTotal();
}
