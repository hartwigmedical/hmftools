package com.hartwig.hmftools.apiclients.diseaseontology.data;

import static java.util.stream.Collectors.toList;

import java.util.List;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Gson.TypeAdapters
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class DiseaseOntologyMetadata {
    @NotNull
    public abstract String id();

    @NotNull
    public abstract String name();

    @Nullable
    public abstract String definition();

    @NotNull
    public abstract List<String> xrefs();

    @NotNull
    public abstract List<List<String>> children();

    @NotNull
    public abstract List<List<String>> parents();

    @NotNull
    @Value.Lazy
    public List<String> parentDoids() {
        return parents().stream().flatMap(List::stream).filter(string -> string.startsWith("DOID:")).collect(toList());
    }

    @NotNull
    @Value.Lazy
    public List<String> childrenDoids() {
        return children().stream().flatMap(List::stream).filter(string -> string.startsWith("DOID:")).collect(toList());
    }
}
