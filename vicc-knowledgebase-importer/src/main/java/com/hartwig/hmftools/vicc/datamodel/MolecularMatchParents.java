package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchParents {

    @NotNull
    public abstract List<String> transcripts();

    @Nullable
    public abstract String type();

    @NotNull
    public abstract String name();

    @Nullable
    public abstract String actionableParent();

}
