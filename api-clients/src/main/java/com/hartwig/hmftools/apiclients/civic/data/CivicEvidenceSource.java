package com.hartwig.hmftools.apiclients.civic.data;

import com.google.gson.annotations.SerializedName;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Gson.TypeAdapters
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicEvidenceSource {
    @Nullable
    public abstract String name();

    public abstract String citation();

    @SerializedName("source_url")
    public abstract String url();
}
