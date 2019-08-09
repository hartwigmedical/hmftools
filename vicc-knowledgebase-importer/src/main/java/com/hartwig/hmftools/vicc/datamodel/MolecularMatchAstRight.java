package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchAstRight {

    @Nullable
    public abstract String operator();

    @Nullable
    public abstract MolecularMatchAstRightRight right();

    @Nullable
    public abstract MolecularMatchAstRightLeft left();

    @Nullable
    public abstract String raw();

    @NotNull
    public abstract String type();

    @Nullable
    public abstract String value();
}
