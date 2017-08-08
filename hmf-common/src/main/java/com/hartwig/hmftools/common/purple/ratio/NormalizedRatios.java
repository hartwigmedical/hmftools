package com.hartwig.hmftools.common.purple.ratio;

import java.util.List;

import com.google.common.collect.Multimap;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface NormalizedRatios {
    @NotNull
    Multimap<String, ReadRatio> normalisedRatios();

    @NotNull
    List<GCMedian> medianReadCount();
}
