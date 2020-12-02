package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTierExplanation {

    @NotNull
    public abstract String tier();

    @NotNull
    public abstract String step();

    @NotNull
    public abstract String message();

    @NotNull
    public abstract String success();
}
