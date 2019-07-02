package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchGRch37Location {

    @NotNull
    public abstract String compositeKey();

    @NotNull
    public abstract String ref();

    @NotNull
    public abstract String stop();

    @NotNull
    public abstract String start();

    @NotNull
    public abstract String chr();

    @NotNull
    public abstract String alt();

    @NotNull
    public abstract String validated();

    @NotNull
    public abstract List<MolecularMatchTranscriptConsequencesGRCH37> transcriptConsequences();

    @NotNull
    public abstract String strand();

}
