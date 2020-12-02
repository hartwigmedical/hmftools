package com.hartwig.hmftools.vicc.datamodel.molecularmatchtrials;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTrialsGeo {

    @NotNull
    public abstract String lat();

    @NotNull
    public abstract String lon();
}
