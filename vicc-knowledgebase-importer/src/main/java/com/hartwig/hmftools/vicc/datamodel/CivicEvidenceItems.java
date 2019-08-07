package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicEvidenceItems {

    @NotNull
    public abstract String status();

    @Nullable
    public abstract String rating();

    @Nullable
    public abstract String drugInteractionType();

    @NotNull
    public abstract String description();

    @NotNull
    public abstract String openChangeCount();

    @NotNull
    public abstract String evidenceType();

    @NotNull
    public abstract List<CivicDrugs> drugs();

    @Nullable
    public abstract String variantOrigin();

    @NotNull
    public abstract CivicDisease disease();

    @NotNull
    public abstract CivicSource source();

    @Nullable
    public abstract String evidenceDirection();

    @NotNull
    public abstract String variantId();

    @Nullable
    public abstract String clinicalSignificance();

    @NotNull
    public abstract String evidenceLevel();

    @NotNull
    public abstract String type();

    @NotNull
    public abstract String id();

    @NotNull
    public abstract String name();




}
