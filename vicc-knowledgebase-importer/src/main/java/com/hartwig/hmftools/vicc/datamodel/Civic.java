package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Civic implements KbSpecificObject {

    @NotNull
    public abstract List<String> variantGroups();

    @NotNull
    public abstract String entrezName();

    @NotNull
    public abstract List<CivicVariantTypes> variantTypes();

    @Nullable
    public abstract String civicActionabilityScore();

    @NotNull
    public abstract List<String> clinvarEntries();

    @NotNull
    public abstract CivicLifecycleActions lifecycleActions();

    @NotNull
    public abstract List<String> variantAliases();

    @Nullable
    public abstract String alleleRegistryId();

    @NotNull
    public abstract String geneId();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract List<CivicEvidenceItems> evidenceItem();

    @NotNull
    public abstract String entrezId();

    @NotNull
    public abstract List<String> assertions();

    @NotNull
    public abstract List<String> hgvs_expressions();

    @NotNull
    public abstract CivicCoordinates coordinates();

    @NotNull
    public abstract String type();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String description();
}
