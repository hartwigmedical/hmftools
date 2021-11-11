package com.hartwig.hmftools.common.lims;

import com.google.gson.annotations.SerializedName;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Gson.TypeAdapters
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
abstract class LimsJsonSubmissionData {

    @NotNull
    @SerializedName("submission")
    public abstract String submission();

    @NotNull
    @SerializedName("project_name")
    public abstract String projectName();

    @NotNull
    @SerializedName("requester_report_contact_email")
    public abstract String reportContactEmail();

    @NotNull
    @SerializedName("requester_report_contact_name")
    public abstract String reportContactName();
}
