package com.hartwig.hmftools.common.flagstat;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Flagstat {

    public abstract long uniqueReadCount();

    public abstract long secondaryCount();

    public abstract long supplementaryCount();

    public abstract double duplicateProportion();

    public abstract double mappedProportion();

    public abstract double pairedInSequencingProportion();

    public abstract double properlyPairedProportion();

    public abstract double withItselfAndMateMappedProportion();

    public abstract double singletonProportion();

}
