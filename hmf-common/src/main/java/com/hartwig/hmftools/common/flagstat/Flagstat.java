package com.hartwig.hmftools.common.flagstat;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Flagstat {

    public abstract double mappedProportion();

    public abstract double duplicateProportion();
}
