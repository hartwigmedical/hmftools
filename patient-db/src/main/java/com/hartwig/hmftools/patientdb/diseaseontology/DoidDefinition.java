package com.hartwig.hmftools.patientdb.diseaseontology;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class DoidDefinition {

    @Nullable
    public abstract String definitionVal();

    @NotNull
    public abstract List<String> definitionXref();
}
