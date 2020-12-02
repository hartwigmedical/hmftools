package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Taxonomy {

    @NotNull
    public abstract String kingdom();

    @NotNull
    public abstract String directParent();

    // Extra S to avoid name-clash with java built-in class.
    @NotNull
    public abstract String classs();

    @Nullable
    public abstract String subClass();

    @NotNull
    public abstract String superClass();
}
