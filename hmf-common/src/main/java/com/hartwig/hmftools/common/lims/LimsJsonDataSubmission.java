package com.hartwig.hmftools.common.lims;

import com.google.gson.annotations.SerializedName;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Gson.TypeAdapters
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
abstract class LimsJsonDataSubmission {

    @NotNull
    @SerializedName("submission")
    public abstract String submission();

    @NotNull
    @SerializedName("contact_email")
    public abstract String contactEmail();

    @NotNull
    @SerializedName("contact_name")
    public abstract String contactName();

}
