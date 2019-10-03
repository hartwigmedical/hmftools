package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchGRch37Location {

    @NotNull
    public abstract String compositeKey();

    @Nullable
    public abstract String ref();

    @Nullable
    public abstract String stop();

    @Nullable
    public abstract String start();

    @Nullable
    public abstract String chr();

    @Nullable
    public abstract String alt();

    @NotNull
    public abstract String validated();

    @NotNull
    public abstract List<MolecularMatchTranscriptConsequencesGRCH37> transcriptConsequences();

    @NotNull
    public abstract String strand();

}
