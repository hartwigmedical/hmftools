package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class EnvironmentalContext {

    @Nullable
    public abstract String term();

    @NotNull
    public abstract String description();

    @Nullable
    public abstract Taxonomy taxonomy();

    @Nullable
    public abstract String source();

    @Nullable
    public abstract String usanStem();

    @Nullable
    public abstract String toxicity();

    @NotNull
    public abstract List<String> approvedCountries();

    @Nullable
    public abstract String id();

}
