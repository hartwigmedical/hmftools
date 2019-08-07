package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTranscriptConsequence {

    @Nullable
    public abstract String aminoAcidChange();

    @NotNull
    public abstract String compositeKey();

    @Nullable
    public abstract String intronNumber();

    @Nullable
    public abstract String exonNumber();

    @NotNull
    public abstract String suppress();

    @Nullable
    public abstract String stop();

    @NotNull
    public abstract String custom();

    @Nullable
    public abstract String start();

    @Nullable
    public abstract String chr();

    @NotNull
    public abstract String strand();

    @NotNull
    public abstract String validated();

    @NotNull
    public abstract String transcript();

    @Nullable
    public abstract String cdna();

    @NotNull
    public abstract String referenceGenome();

    @Nullable
    public abstract String ref();

    @Nullable
    public abstract String alt();
}
