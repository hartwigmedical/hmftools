package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchAstRight {

    @NotNull
    public abstract String type();

    @Nullable
    public abstract String raw();

    @Nullable
    public abstract String value();

    @Nullable
    public abstract String operator();

    @Nullable
    public abstract MolecularMatchAstRightLeft left();

    @Nullable
    public abstract MolecularMatchAstRightRight right();
}
