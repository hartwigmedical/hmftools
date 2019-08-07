package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchAstLeft {

    @Nullable
    public abstract String operator();

    @Nullable
    public abstract String raw();

    @NotNull
    public abstract String type();

    @NotNull
    public abstract String value();
}
