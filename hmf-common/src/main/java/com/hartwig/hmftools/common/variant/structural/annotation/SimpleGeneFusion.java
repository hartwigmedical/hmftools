package com.hartwig.hmftools.common.variant.structural.annotation;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class SimpleGeneFusion {

    @NotNull
    public abstract String fiveGene();

    @NotNull
    public abstract String threeGene();

}
