package com.hartwig.hmftools.common.doid;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class Doid {

    @NotNull
    public abstract String doid();

    @NotNull
    public abstract String id();

    @Nullable
    public abstract String doidTerm();

    @Nullable
    public abstract String type();

    @Nullable
    public abstract DoidMetadata doidMetadata();
}
