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
abstract class LimsJsonSampleData {

    @NotNull
    @SerializedName("sample_name")
    public abstract String sampleId();

    @NotNull
    @SerializedName("patient")
    public abstract String patientId();

    // ref sample id could be null in lims if it is a reference sample
    @Nullable
    @SerializedName("ref_sample_id")
    public abstract String refBarcodeId();

    @NotNull
    @SerializedName("sample_id")
    public abstract String tumorBarcodeId();

    // Patient number is only used for CORE project
    @Nullable
    @SerializedName("hospital_patient_id")
    public abstract String hospitalPatientId();

    @NotNull
    @SerializedName("arrival_date")
    public abstract String arrivalDateString();

    @NotNull
    @SerializedName("dna_conc")
    public abstract String dnaConcentration();

    // Sampling date is only known for CPCT/DRUP/WIDE/CORE tumor biopsies.
    @Nullable
    @SerializedName("sampling_date")
    public abstract String samplingDateString();

    // Tumor biopsies analyzed in research context do not have a proper tumor percentage filled in.
    @Nullable
    @SerializedName("tumor_perc")
    public abstract String tumorPercentageString();

    @NotNull
    @SerializedName("ptum")
    public abstract String primaryTumor();

    @NotNull
    @SerializedName("lab_sop_versions")
    abstract String labSopVersions();

    // Lab remarks is an optional field in LIMS
    @Nullable
    @SerializedName("lab_remarks")
    public abstract String labRemarks();

    @NotNull
    @SerializedName("label")
    public abstract String labelSample();

    @NotNull
    @SerializedName("project_name")
    public abstract String projectName();

    @NotNull
    @SerializedName("submission")
    public abstract String submission();

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
