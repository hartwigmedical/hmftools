package com.hartwig.hmftools.patientdb.data;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class CuratedTreatment {
    @NotNull
    public abstract String name();

    @NotNull
    public abstract String type();
}
