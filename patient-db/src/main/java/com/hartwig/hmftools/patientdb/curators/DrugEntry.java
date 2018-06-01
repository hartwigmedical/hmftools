package com.hartwig.hmftools.patientdb.curators;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
abstract class DrugEntry {

    @NotNull
    public abstract String canonicalName();

    @NotNull
    public abstract List<String> synonyms();

    @NotNull
    public abstract String type();

    @NotNull
    public abstract String treatmentMechanism();
}
