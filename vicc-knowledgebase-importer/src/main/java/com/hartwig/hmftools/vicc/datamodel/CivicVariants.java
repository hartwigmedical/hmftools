package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicVariants {

    @NotNull
    public abstract String entrez_name();

    @NotNull
    public abstract String variant_types();

    @NotNull
    public abstract String description();

    @NotNull
    public abstract String civic_actionability_score();

    @NotNull
    public abstract String gene_id();

    @NotNull
    public abstract String entrez_id();

    @NotNull
    public abstract String coordinates();

    @NotNull
    public abstract String type();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String name();
}
