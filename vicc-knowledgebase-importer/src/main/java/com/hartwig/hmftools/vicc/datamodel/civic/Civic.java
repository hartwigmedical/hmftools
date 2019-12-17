package com.hartwig.hmftools.vicc.datamodel.civic;

import java.util.List;

import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Civic implements KbSpecificObject {

    @NotNull
    public abstract String entrezId();

    @NotNull
    public abstract String entrezName();

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String type();

    @NotNull
    public abstract CivicCoordinates coordinates();

    @NotNull
    public abstract List<CivicSource> sources();

    @NotNull
    public abstract List<String> variantAliases();

    @NotNull
    public abstract List<CivicVariantGroup> variantGroups();

    @NotNull
    public abstract List<CivicVariantType> variantTypes();

    @NotNull
    public abstract List<String> hgvsExpressions();

    @NotNull
    public abstract CivicEvidenceItem evidenceItem();

    @NotNull
    public abstract List<String> assertions();

    @Nullable
    public abstract String civicActionabilityScore();

    @NotNull
    public abstract List<String> clinVarEntries();

    @Nullable
    public abstract String alleleRegistryId();

    @Nullable
    public abstract CivicProvisionalValue provisionalValue();

    @NotNull
    public abstract CivicLifecycleActions lifecycleActions();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String geneId();

    @NotNull
    public abstract String description();
}
