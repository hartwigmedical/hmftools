package com.hartwig.hmftools.protect.linx;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ViralInsertion {

    @NotNull
    public abstract String virus();

    public abstract int viralInsertionCount();
}
