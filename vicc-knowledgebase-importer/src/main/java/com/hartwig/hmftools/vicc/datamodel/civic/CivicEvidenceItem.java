package com.hartwig.hmftools.vicc.datamodel.civic;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicEvidenceItem {

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String type();

    @NotNull
    public abstract String status();

    @Nullable
    public abstract String rating();

    @NotNull
    public abstract String evidenceType();

    @NotNull
    public abstract String evidenceLevel();

    @Nullable
    public abstract String evidenceDirection();

    @Nullable
    public abstract String drugInteractionType();

    @NotNull
    public abstract List<CivicDrug> drugs();

    @NotNull
    public abstract CivicDisease disease();

    @Nullable
    public abstract String variantOrigin();

    @NotNull
    public abstract CivicSource source();

    @Nullable
    public abstract String clinicalSignificance();

    @NotNull
    public abstract String openChangeCount();

    @NotNull
    public abstract String description();

    @Nullable
    public abstract String variantId();

    @NotNull
    public abstract String id();

}
