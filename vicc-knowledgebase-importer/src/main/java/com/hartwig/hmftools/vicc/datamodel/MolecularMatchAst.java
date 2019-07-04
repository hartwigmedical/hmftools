package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchAst {

    @NotNull
    public abstract String operator();

    @NotNull
    public abstract MolecularMatchAstRight right();

    @NotNull
    public abstract String type();

    @NotNull
    public abstract MolecularMatchAstLeft left();
}
