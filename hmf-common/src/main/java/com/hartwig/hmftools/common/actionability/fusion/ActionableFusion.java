package com.hartwig.hmftools.common.actionability.fusion;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class ActionableFusion {

    @NotNull
    public abstract String fiveGene();

    @NotNull
    public abstract String threeGene();

    @NotNull
    public abstract String source();

    @NotNull
    public abstract String reference();

    @NotNull
    public abstract String drug();

    @NotNull
    public abstract String drugsType();

    @NotNull
    public abstract String cancerType();

    @NotNull
    public abstract String level();

    @NotNull
    public abstract String response();
}
