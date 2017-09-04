package com.hartwig.hmftools.apiclients.civic.data;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Gson.TypeAdapters
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicDrug {
    public abstract String id();

    public abstract String name();

    @Override
    public String toString() {
        return name();
    }
}
