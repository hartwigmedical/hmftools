package com.hartwig.hmftools.patientdb.diseaseontology;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DoidMetadata {

    @NotNull
    public abstract DoidDefinition doidDefinition();

    @NotNull
    public abstract List<String> subset();

    @NotNull
    public abstract List<DoidXref> xref();

    @NotNull
    public abstract List<DoidSynonyms> synonyms();

    @NotNull
    public abstract List<DoidBasicPropertyValues> basicPropertyValues();
}
