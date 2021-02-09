package com.hartwig.hmftools.ckb.datamodelinterpretation.gene;

import java.util.Date;
import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Gene {

    public abstract int id();

    @NotNull
    public abstract String geneSymbol();

    @NotNull
    public abstract List<String> terms();

    @Nullable
    public abstract String entrezId();

    @NotNull
    public abstract List<String> synonyms();

    @Nullable
    public abstract String chromosome();

    @Nullable
    public abstract String mapLocation();

//    @NotNull
//    public abstract List<GeneDescription> geneDescriptions();
//
//    @NotNull
//    public abstract String canonicalTranscript();
//
//    @NotNull
//    public abstract String geneRole();
//
//    @NotNull
//    public abstract Date createDate();
//
//    @NotNull
//    public abstract Date updateDate();
//
//    @NotNull
//    public abstract List<GeneClinicalTrial> geneClinicalTrials();
//
//    @NotNull
//    public abstract List<GeneEvidence> geneEvidences();
}