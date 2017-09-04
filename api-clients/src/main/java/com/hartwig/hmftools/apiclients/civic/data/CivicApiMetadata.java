package com.hartwig.hmftools.apiclients.civic.data;

import com.google.gson.annotations.SerializedName;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Gson.TypeAdapters
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicApiMetadata {

    @SerializedName("current_page")
    public abstract int currentPage();

    @SerializedName("per_page")
    public abstract int perPage();

    @SerializedName("total_pages")
    public abstract int totalPages();

    @SerializedName("total_count")
    public abstract int totalElements();
}
