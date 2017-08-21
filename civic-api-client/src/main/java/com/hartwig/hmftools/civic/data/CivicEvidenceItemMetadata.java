package com.hartwig.hmftools.civic.data;

import com.google.gson.annotations.SerializedName;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Gson.TypeAdapters
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class CivicEvidenceItemMetadata {
    @SerializedName("accepted_count")
    public abstract int accepted();

    @SerializedName("rejected_count")
    public abstract int rejected();

    @SerializedName("submitted_count")
    public abstract int submitted();
}
