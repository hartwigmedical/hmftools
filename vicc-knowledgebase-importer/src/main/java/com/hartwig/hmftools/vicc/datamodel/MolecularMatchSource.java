package com.hartwig.hmftools.vicc.datamodel;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class MolecularMatchSource {

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String suppress();

    @NotNull
    public abstract String pubId();

    @Nullable
    public abstract String subType();

    @NotNull
    public abstract String valid();

    @NotNull
    public abstract String link();

    @NotNull
    public abstract String year();

    @Nullable
    public abstract String trialId();

    @NotNull
    public abstract String type();

    @NotNull
    public abstract String id();

    @Nullable
    public abstract String institution();

    @Nullable
    public abstract String trialPhase();

    @Nullable
    public abstract String functionalConsequence();

    @Nullable
    public abstract String trustRating();
}
