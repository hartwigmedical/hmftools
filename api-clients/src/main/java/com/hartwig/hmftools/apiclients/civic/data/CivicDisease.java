package com.hartwig.hmftools.apiclients.civic.data;

import com.google.gson.annotations.SerializedName;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Gson.TypeAdapters
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicDisease implements Comparable<CivicDisease> {
    public abstract int id();

    public abstract String name();

    @SerializedName("display_name")
    public abstract String displayName();

    @Nullable
    @SerializedName("doid")
    public abstract String ontologyId();

    @NotNull
    @Value.Derived
    public String doidString() {
        return "DOID:" + ontologyId();
    }

    @Override
    public int compareTo(@NotNull final CivicDisease other) {
        return Integer.compare(id(), other.id());
    }
}
