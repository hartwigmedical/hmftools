package com.hartwig.hmftools.orange.algo.wildtype;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class WildTypeGene {

    @NotNull
    public abstract String gene();
}