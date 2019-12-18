package com.hartwig.hmftools.vicc.datamodel.molecularmatch;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchMutation {

    @NotNull
    public abstract String geneSymbol();

    @NotNull
    public abstract String name();

    @Nullable
    public abstract String transcriptRecognized();

    @Nullable
    public abstract String transcript();

    @Nullable
    public abstract String longestTranscript();

    @Nullable
    public abstract String uniprotTranscript();

    @Nullable
    public abstract List<MolecularMatchTranscriptConsequence> transcriptConsequences();

    @Nullable
    public abstract List<MolecularMatchParent> parents();

    @Nullable
    public abstract List<MolecularMatchWGSAData> wgsaData();

    @Nullable
    public abstract List<MolecularMatchWGSaMap2> wgsaMaps();

    @Nullable
    public abstract MolecularMatchExonsInfo exonsInfo();

    @Nullable
    public abstract List<MolecularMatchFusionData> fusionData();

    @NotNull
    public abstract List<String> mutationTypes();

    @NotNull
    public abstract List<String> sources();

    @NotNull
    public abstract List<String> synonyms();

    @NotNull
    public abstract List<MolecularMatchGRch37Location2> grch37Locations();

    @NotNull
    public abstract List<String> pathology();

    @NotNull
    public abstract List<String> cDNA();

    @NotNull
    public abstract String description();

    @NotNull
    public abstract String src();

    @NotNull
    public abstract String id();
}
