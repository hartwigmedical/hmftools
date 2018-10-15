package com.hartwig.hmftools.common.actionability.fusion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class ActionabilityFusionPairs {

    @NotNull
    public abstract String fiveGene();

    @NotNull
    public abstract String threeGene();

    @NotNull
    public abstract String source();

    @NotNull
    public abstract String reference();

    @NotNull
    public abstract String drugsName();

    @NotNull
    public abstract String drugsType();

    @NotNull
    public abstract String cancerType();

    @NotNull
    public abstract String hmfLevel();

    @NotNull
    public abstract String hmfResponse();
}
