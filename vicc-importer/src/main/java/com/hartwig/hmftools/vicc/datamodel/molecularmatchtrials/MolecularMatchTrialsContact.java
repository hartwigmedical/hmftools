package com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTrialsContact {

    @Nullable
    public abstract String name();

    @Nullable
    public abstract String email();

    @Nullable
    public abstract String phone();
}
