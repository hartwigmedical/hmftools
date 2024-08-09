package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.common.purple.GermlineStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Segment
{
    public abstract int bafCount();

    public abstract double observedTumorRatio();

    @NotNull
    public abstract GermlineStatus germlineStatus();
}