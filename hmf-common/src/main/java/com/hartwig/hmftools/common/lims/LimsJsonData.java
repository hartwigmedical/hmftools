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
abstract class LimsJsonData {

    @NotNull
    @SerializedName("sample_name")
    public abstract String sampleId();

    @NotNull
    @SerializedName("arrival_date")
    public abstract String arrivalDateString();

    @NotNull
    @SerializedName("dna_conc")
    public abstract String dnaConcentration();

    // KODU: Sampling date is only known for CPCT/DRUP tumor biopsies.
    @Nullable
    @SerializedName("sampling_date")
    public abstract String samplingDateString();

    // KODU: Tumor biopsies analyzed in research context do not have a proper tumor percentage filled in.
    // LISC: Tumor biopsies could be NA as tumor percentage
    @Nullable
    @SerializedName("tumor_perc")
    public abstract String tumorPercentageString();

    @NotNull
    @SerializedName("ptum")
    public abstract String primaryTumor();

    @NotNull
    @SerializedName("lab_sop_versions")
    abstract String labSopVersions();

    @Nullable
    @SerializedName("lab_remarks")
    public abstract String labRemarks();

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
