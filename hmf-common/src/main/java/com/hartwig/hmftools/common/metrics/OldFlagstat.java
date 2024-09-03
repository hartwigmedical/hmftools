package com.hartwig.hmftools.common.metrics;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Deprecated
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class OldFlagstat
{
    public abstract long uniqueReadCount();

    public abstract long secondaryCount();

    public abstract long supplementaryCount();

    public abstract double duplicateProportion();

    public abstract double mappedProportion();

    public abstract double pairedInSequencingProportion();

    public abstract double properlyPairedProportion();

    public abstract double withItselfAndMateMappedProportion();

    public abstract double singletonProportion();

    public static final double MIN_MAPPED_PROPORTION = 0.95;
    public static boolean pass(@NotNull OldFlagstat flagstat)
    {
        return pass(flagstat.mappedProportion());
    }
    public static boolean pass(double mappedProportion)
    {
        return mappedProportion >= MIN_MAPPED_PROPORTION;
    }

}
