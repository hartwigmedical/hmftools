package com.hartwig.hmftools.patientdb.data;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class CuratedTreatment implements Comparable<CuratedTreatment> {
    @NotNull
    public abstract String name();

    @NotNull
    public abstract String type();

    @NotNull
    public abstract String searchTerm();

    @Override
    public int compareTo(@NotNull final CuratedTreatment other) {
        return name().compareTo(other.name());
    }
}
