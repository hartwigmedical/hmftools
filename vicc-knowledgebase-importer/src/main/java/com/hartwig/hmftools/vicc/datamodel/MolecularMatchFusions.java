package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchFusions {

    @NotNull
    public abstract String referenceGenome();

    @NotNull
    public abstract String LBPWREP();

    @NotNull
    public abstract String RBPWREP();

    @NotNull
    public abstract String exonNumber();

    @NotNull
    public abstract String chr();

    @NotNull
    public abstract String RBPWLEP();

    @NotNull
    public abstract String intronNumber();

    @NotNull
    public abstract String LBPWLEP();
}
