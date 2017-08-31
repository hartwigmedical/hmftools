package com.hartwig.hmftools.common.lims;

import com.google.gson.annotations.SerializedName;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Gson.TypeAdapters
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class LimsJsonData {
    @NotNull
    @SerializedName("sample_source")
    public abstract String sampleSource();

    @NotNull
    public abstract String patient();

    @NotNull
    @SerializedName("sample_name")
    public abstract String sampleName();

    @NotNull
    @SerializedName("sample_barcode")
    public abstract String sampleBarcode();

    @NotNull
    @SerializedName("arrival_date")
    public abstract String arrivalDateString();

    @NotNull
    @SerializedName("sampling_date")
    public abstract String samplingDateString();

    @NotNull
    @SerializedName("tumor_perc")
    public abstract String tumorPercentage();
}
