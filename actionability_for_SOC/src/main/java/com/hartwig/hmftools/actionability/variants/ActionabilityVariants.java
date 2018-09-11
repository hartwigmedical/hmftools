package com.hartwig.hmftools.actionability.variants;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })

public abstract class ActionabilityVariants {

    @NotNull
    abstract String gene();

    @NotNull
    abstract String chromosome();

    @NotNull
    abstract String chromosmePosition();


}
