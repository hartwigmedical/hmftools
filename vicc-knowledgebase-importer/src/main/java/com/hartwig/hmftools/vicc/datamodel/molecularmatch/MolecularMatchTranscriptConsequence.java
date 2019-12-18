package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchTranscriptConsequence {

    @Nullable
    public abstract String chr();

    @Nullable
    public abstract String start();

    @Nullable
    public abstract String stop();

    @Nullable
    public abstract String ref();

    @Nullable
    public abstract String alt();

    @NotNull
    public abstract String referenceGenome();

    @NotNull
    public abstract String transcript();

    @NotNull
    public abstract String strand();

    @Nullable
    public abstract String cdna();

    @Nullable
    public abstract String aminoAcidChange();

    @Nullable
    public abstract String intronNumber();

    @NotNull
    public abstract List<String> exonNumbers();

    @NotNull
    public abstract String suppress();

    @NotNull
    public abstract String custom();

    @NotNull
    public abstract String validated();

    @NotNull
    public abstract String compositeKey();
}
