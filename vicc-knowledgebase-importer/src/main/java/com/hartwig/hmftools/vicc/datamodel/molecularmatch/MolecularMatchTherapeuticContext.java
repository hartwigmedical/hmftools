package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTherapeuticContext {

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String facet();

    @NotNull
    public abstract String suppress();

    @Nullable
    public abstract String valid();
}
