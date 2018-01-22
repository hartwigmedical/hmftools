package com.hartwig.hmftools.common.cosmic.genes;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class CosmicGeneData {

    @NotNull
    public abstract String description();

    @NotNull
    public abstract String entrezId();

    @NotNull
    public abstract String genomeLocation();

    @NotNull
    public abstract String chromosomeBand();

    @NotNull
    public abstract String somatic();

    @NotNull
    public abstract String germline();

    @NotNull
    public abstract String somaticTumorTypes();

    @NotNull
    public abstract String germlineTumorTypes();

    @NotNull
    public abstract String cancerSyndrome();

    @NotNull
    public abstract String tissueType();

    @NotNull
    public abstract String molecularGenetics();

    @NotNull
    public abstract String role();

    @NotNull
    public abstract String mutationTypes();

    @NotNull
    public abstract String translocationPartner();

    @NotNull
    public abstract String otherGermlineMut();

    @NotNull
    public abstract String otherSyndrome();
}
