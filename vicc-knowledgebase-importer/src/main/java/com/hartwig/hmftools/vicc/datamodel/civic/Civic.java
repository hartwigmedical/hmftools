package com.hartwig.hmftools.vicc.datamodel.civic;

import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Civic implements KbSpecificObject {

    @Nullable
    public abstract List<CivicVariantGroup> variantGroups();

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

    @Nullable
    public abstract CivicDescription provisional_values();

    @NotNull
    public abstract String geneId();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract List<CivicEvidenceItems> evidenceItem();

    @Nullable
    public abstract List<CivicSource> sources();

    @NotNull
    public abstract String entrezId();

    @NotNull
    public abstract List<String> assertions();

    @NotNull
    public abstract List<String> hgvs_expressions();

    @NotNull
    public abstract CivicError errors();

    @NotNull
    public abstract CivicCoordinates coordinates();

    @NotNull
    public abstract String type();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String description();
}
