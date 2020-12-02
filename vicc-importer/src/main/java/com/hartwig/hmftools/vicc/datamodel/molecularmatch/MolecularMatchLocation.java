package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchLocation {

    @NotNull
    public abstract String chr();

    @NotNull
    public abstract String start();

    @NotNull
    public abstract String stop();

    @Nullable
    public abstract String ref();

    @Nullable
    public abstract String alt();

    @Nullable
    public abstract String cdna();

    @Nullable
    public abstract String aminoAcidChange();

    @Nullable
    public abstract String referenceGenome();

    @Nullable
    public abstract String strand();

    @Nullable
    public abstract String intronNumber();

    @NotNull
    public abstract List<String> exonNumbers();
}
