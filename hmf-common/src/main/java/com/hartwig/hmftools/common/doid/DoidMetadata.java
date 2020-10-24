package com.hartwig.hmftools.common.doid;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DoidMetadata {

    @Nullable
    public abstract DoidDefinition doidDefinition();

    @NotNull
    public abstract List<String> subsets();

    @NotNull
    public abstract List<DoidXref> xrefs();

    @Nullable
    public abstract List<DoidSynonym> synonyms();

    @Nullable
    public abstract List<DoidBasicPropertyValue> basicPropertyValues();
}
