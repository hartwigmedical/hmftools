package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Association {

    @NotNull
    public abstract List<String> variantNames();

    @NotNull
    public abstract Evidence evidence();

    @Nullable
    public abstract String evidenceLevel();

    @Nullable
    public abstract String evidenceLabel();

    @Nullable
    public abstract String responseType();

    @Nullable
    public abstract String drugLabels();

    @Nullable
    public abstract String sourceLink();

    @NotNull
    public abstract List<String> publicationUrls();

    @Nullable
    public abstract Phenotype phenotype();

    @NotNull
    public abstract String description();

    @Nullable
    public abstract List<EnvironmentalContext> environmentalContexts();

    @Nullable
    public abstract String oncogenic();
}
