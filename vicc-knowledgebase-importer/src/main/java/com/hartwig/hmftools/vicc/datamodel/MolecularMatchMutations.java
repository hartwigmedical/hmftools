package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchMutations {

    @Nullable
    public abstract List<MolecularMatchTranscriptConsequence> transcriptConsequence();

    @NotNull
    public abstract String description();

    @NotNull
    public abstract List<String> mutationType();

    @NotNull
    public abstract String src();

    @NotNull
    public abstract List<String> sources();

    @NotNull
    public abstract List<String> synonyms();

    @NotNull
    public abstract List<MolecularMatchGRch37Location> gRch37Location();

    @Nullable
    public abstract String uniprotTranscript();

    @NotNull
    public abstract String geneSymbol();

    @NotNull
    public abstract List<String> pathology();

    @Nullable
    public abstract String transcript();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract List<String> cDNA();

    @NotNull
    public abstract String name();
}
