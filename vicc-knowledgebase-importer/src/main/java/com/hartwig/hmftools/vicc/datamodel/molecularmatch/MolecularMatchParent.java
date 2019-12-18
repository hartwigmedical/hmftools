package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchParent {

    @NotNull
    public abstract String name();

    @Nullable
    public abstract String type();

    @Nullable
    public abstract String actionableParent();

    @NotNull
    public abstract List<String> transcripts();

}
