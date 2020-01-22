package com.hartwig.hmftools.protect.conclusion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class TumorLocationConclusion {

    @NotNull
    public abstract String primaryTumorLocation();

    @NotNull
    public abstract String cancerSubType();

    @NotNull
    public abstract String tumorLocationConclusion();
}
