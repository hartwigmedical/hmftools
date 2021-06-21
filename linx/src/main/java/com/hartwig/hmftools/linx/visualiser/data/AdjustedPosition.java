package com.hartwig.hmftools.linx.visualiser.data;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class AdjustedPosition implements GenomePosition
{
    public abstract long unadjustedPosition();

    public abstract int svId();
}
