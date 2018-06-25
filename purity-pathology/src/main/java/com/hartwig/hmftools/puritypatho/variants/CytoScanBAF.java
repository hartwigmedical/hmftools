package com.hartwig.hmftools.puritypatho.variants;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class CytoScanBAF {
    @NotNull
    public abstract int chromosomeNumber();

    @NotNull
    public abstract int postionStart();

    @NotNull
    public abstract int positionEnd();
}
