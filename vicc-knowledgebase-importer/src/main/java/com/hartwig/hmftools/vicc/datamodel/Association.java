package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Association {

    @NotNull
    public abstract String variantName();

    @NotNull
    public abstract List<Evidence> evidence();

    @NotNull
    public abstract String evidenceLevel();

    @NotNull
    public abstract String evidenceLabel();

    @NotNull
    public abstract String responseType();

    @NotNull
    public abstract String drugLabels();

    @NotNull
    public abstract String sourceLink();

    @NotNull
    public abstract List<String> publicationUrls();

    @NotNull
    public abstract Phenotype phenotype();

    @NotNull
    public abstract String description();

    @NotNull
    public abstract List<EnvironmentalContext> environmentalContexts();
}
