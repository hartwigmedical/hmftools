package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class ViccEntry {

    @NotNull
    public abstract String source();

    @NotNull
    public abstract List<String> genes();

    @NotNull
    public abstract List<GeneIdentifier> geneIdentifiers();

    // Field does not exist for SAGE records
    @Nullable
    public abstract List<String> featureNames();

    @NotNull
    public abstract List<Feature> features();

    @NotNull
    public abstract Association association();

    @NotNull
    public abstract List<String> tags();

    @NotNull
    public abstract List<String> devTags();

    @NotNull
    public abstract Cgi cgi();

    @NotNull
    public abstract BRCA brca();

    @NotNull
    public abstract Sage sage();

    @NotNull
    public abstract Pmkb pmkb();

}

