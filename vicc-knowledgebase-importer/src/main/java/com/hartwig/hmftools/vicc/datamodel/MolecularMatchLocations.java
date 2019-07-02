package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchLocations {

    @NotNull
    public abstract String aminoAcidChange();

    @NotNull
    public abstract String intronNumber();

    @NotNull
    public abstract String exonNumber();

    @NotNull
    public abstract String stop();

    @NotNull
    public abstract String start();

    @NotNull
    public abstract String chr();

    @NotNull
    public abstract String strand();

    @NotNull
    public abstract String alt();

    @NotNull
    public abstract String referenceGenome();

    @NotNull
    public abstract String ref();

    @NotNull
    public abstract String cdna();
}
