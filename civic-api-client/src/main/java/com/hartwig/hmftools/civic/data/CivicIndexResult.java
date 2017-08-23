package com.hartwig.hmftools.civic.data;

import java.util.List;

import com.google.gson.annotations.SerializedName;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Gson.TypeAdapters
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicIndexResult<T> {

    @SerializedName("_meta")
    public abstract CivicApiMetadata meta();

    public abstract List<T> records();
}
