package com.hartwig.hmftools.vicc.datamodel.civic;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicVariant {

    @NotNull
    public abstract String entrezId();

    @NotNull
    public abstract String entrezName();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String type();

    @NotNull
    public abstract List<CivicVariantType> variantTypes();

    @Nullable
    public abstract String civicActionabilityScore();

    @Nullable
    public abstract CivicCoordinates coordinates();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String geneId();

    @NotNull
    public abstract String description();
}
