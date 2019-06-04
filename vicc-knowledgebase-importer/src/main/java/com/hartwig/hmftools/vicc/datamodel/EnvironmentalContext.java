package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class EnvironmentalContext {

    @NotNull
    public abstract String term();

    @NotNull
    public abstract String description();

    @NotNull
    public abstract Taxonomy taxonomy();

    @NotNull
    public abstract String source();

    @NotNull
    public abstract String usanStem();

    @NotNull
    public abstract List<String> approvedCountries();

    @NotNull
    public abstract String id();

}
