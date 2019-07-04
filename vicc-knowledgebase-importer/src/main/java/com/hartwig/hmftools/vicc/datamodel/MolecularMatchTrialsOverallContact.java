package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTrialsOverallContact {

    @NotNull
    public abstract String phone();

    @NotNull
    public abstract String last_name();

    @NotNull
    public abstract String email();

    @NotNull
    public abstract String affiliation();
}
