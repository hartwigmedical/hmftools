package com.hartwig.hmftools.linx.visualiser.data;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Connector implements GenomePosition
{
    public abstract int frame();

    public abstract int clusterId();

    public abstract int chainId();

    public abstract int track();

    public abstract double ploidy();
}
