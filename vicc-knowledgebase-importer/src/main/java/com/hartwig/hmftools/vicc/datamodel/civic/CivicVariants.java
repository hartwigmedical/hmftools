package com.hartwig.hmftools.vicc.datamodel.civic;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicVariants {

    @NotNull
    public abstract String entrez_name();

    @NotNull
    public abstract List<CivicVariantTypes> variant_types();

    @NotNull
    public abstract String description();

    @Nullable
    public abstract String civic_actionability_score();

    @NotNull
    public abstract String gene_id();

    @NotNull
    public abstract String entrez_id();

    @Nullable
    public abstract CivicCoordinates coordinates();

    @NotNull
    public abstract String type();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String name();
}
