package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchAst {

    @Nullable
    public abstract String raw();

    @Nullable
    public abstract String value();

    @NotNull
    public abstract String operator();

    @Nullable
    public abstract MolecularMatchAstRight right();

    @NotNull
    public abstract String type();

    @Nullable
    public abstract MolecularMatchAstLeft left();
}
