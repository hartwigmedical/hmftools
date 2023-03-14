package com.hartwig.hmftools.datamodel.cuppa;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CuppaPrediction {

    @NotNull
    public abstract String cancerType();

    public abstract double likelihood();
}
