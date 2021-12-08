package com.hartwig.hmftools.serve.sources.actin;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActinTrial {

    @NotNull
    public abstract String trialId();

    @NotNull
    public abstract String cohortId();

    @NotNull
    public abstract String rule();

    @NotNull
    public abstract String parameters();
}
