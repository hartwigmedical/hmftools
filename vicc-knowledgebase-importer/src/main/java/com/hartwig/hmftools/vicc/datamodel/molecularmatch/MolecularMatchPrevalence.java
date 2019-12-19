package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchPrevalence {

    @NotNull
    public abstract String studyId();

    @NotNull
    public abstract String count();

    @NotNull
    public abstract String samples();

    @NotNull
    public abstract String percent();

    @Nullable
    public abstract String molecular();

    @Nullable
    public abstract String condition();
}
