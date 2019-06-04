package com.hartwig.hmftools.vicc.datamodel;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class Feature {

    @NotNull
    public abstract String name();

    @NotNull
    public abstract String biomarkerType();

    @NotNull
    public abstract String referenceName();

    @NotNull
    public abstract String chromosome();

    @NotNull
    public abstract String start();

    @NotNull
    public abstract String end();

    @NotNull
    public abstract String ref();

    @NotNull
    public abstract String alt();

    @NotNull
    public abstract List<String> provenance();

    @NotNull
    public abstract String provenanceRule();

    @NotNull
    public abstract String geneSymbol();

    @NotNull
    public abstract List<String> synonyms();

    @NotNull
    public abstract String entrezId();

    @NotNull
    public abstract SequenceOntology sequenceOntology();

    @NotNull
    public abstract List<String> links();

    @NotNull
    public abstract String description();
}
