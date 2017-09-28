package com.hartwig.hmftools.common.lims;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

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

    @NotNull
    @SerializedName("lab_sop_versions")
    protected abstract String labSopVersions();

    @NotNull
    @Value.Derived
    public String labProcedures() {
        final Pattern pattern = Pattern.compile("PREP(\\d+)V(\\d+)-QC(\\d+)V(\\d+)-SEQ(\\d+)V(\\d+)");
        final Matcher matcher = pattern.matcher(labSopVersions());
        if (matcher.find()) {
            return labSopVersions();
        } else {
            return "N/A";
        }
    }
}
