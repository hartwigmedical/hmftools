package com.hartwig.hmftools.vicc.datamodel.pmkb;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class PmkbVariant {

    @Nullable
    public abstract String name();

    @Nullable
    public abstract String coordinates();

    @Nullable
    public abstract String chromosome();

    @Nullable
    public abstract String cytoband();

    @NotNull
    public abstract PmkbGene gene();

    @Nullable
    public abstract String transcript();

    @Nullable
    public abstract String effect();

    @Nullable
    public abstract String codons();

    @Nullable
    public abstract String exons();

    @Nullable
    public abstract String dnaChange();

    @Nullable
    public abstract String aminoAcidChange();

    @Nullable
    public abstract String germline();

    @Nullable
    public abstract String partnerGene();

    @Nullable
    public abstract String cnvType();

    @Nullable
    public abstract String chromosomeBasedCnv();

    @Nullable
    public abstract String variantType();

    @Nullable
    public abstract String cosmic();

    @Nullable
    public abstract String description();

    @Nullable
    public abstract String descriptionType();

    @Nullable
    public abstract String notes();

    @Nullable
    public abstract String id();
}
