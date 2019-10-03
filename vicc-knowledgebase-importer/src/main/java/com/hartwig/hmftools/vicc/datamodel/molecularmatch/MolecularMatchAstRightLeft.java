package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchAstRightLeft {
    @Nullable
    public abstract String raw();

    @NotNull
    public abstract String type();

    @Nullable
    public abstract String value();

    @Nullable
    public abstract MolecularMatchAstRightLeftRight right();

    @Nullable
    public abstract MolecularMatchAstLeft left();
}
