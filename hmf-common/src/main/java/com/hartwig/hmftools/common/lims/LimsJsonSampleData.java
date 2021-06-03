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

    // ref_sample_id could be null in lims if the sample itself is a reference sample
    @Nullable
    @SerializedName("ref_sample_id")
    public abstract String refBarcode();

    @NotNull
    @SerializedName("sample_id")
    public abstract String tumorBarcode();

    @NotNull
    @SerializedName("arrival_date")
    public abstract String arrivalDate();

    // Sampling date is only known for CPCT/DRUP/WIDE/CORE tumor biopsies.
    @Nullable
    @SerializedName("sampling_date")
    public abstract String samplingDate();

    @SerializedName("shallowseq")
    public abstract boolean shallowSeq();

    // Tumor biopsies analyzed in research context do not have a pathology tumor percentage. Also WIDE samples are no longer sent to PA.
    @Nullable
    @SerializedName("tumor_perc")
    public abstract String pathologyTumorPercentage();

    @NotNull
    @SerializedName("conc")
    public abstract String dnaConcentration();

    @NotNull
    @SerializedName("ptum")
    public abstract String primaryTumor();

    // Biopsy Location is only known for COREDB samples
    @Nullable
    @SerializedName("biopsy_site")
    public abstract String biopsySite();

    @NotNull
    @SerializedName("submission")
    public abstract String submission();

    // Patient number is only used for CORE project
    @Nullable
    @SerializedName("hospital_patient_id")
    public abstract String hospitalPatientId();

    // Pathology sample id is always used for WIDE project and is optional for CORE projects
    @Nullable
    @SerializedName("hospital_pa_sample_id")
    public abstract String hospitalPathologySampleId();

    @SerializedName("report_germline")
    public abstract boolean reportGermlineVariants();

    // Choice regarding reporting of germline findings is only used in WIDE projects
    @Nullable
    @SerializedName("report_germline_level")
    public abstract String germlineReportingLevel();

    @NotNull
    @SerializedName("lab_sop_versions")
    abstract String labSopVersions();

    @SerializedName("report_viral")
    public abstract boolean reportViralPresence();

    @Nullable
    @SerializedName("cohort")
    public abstract String cohort();

    @Nullable
    @SerializedName("analysis_type")
    public abstract String analysisType();

    @NotNull
    @Value.Derived
    public String labProcedures() {
        Pattern pattern = Pattern.compile("PREP(\\d+)V(\\d+)-QC(\\d+)V(\\d+)-SEQ(\\d+)V(\\d+)");
        Matcher matcher = pattern.matcher(labSopVersions());
        if (matcher.find()) {
            return labSopVersions();
        } else {
            return Lims.NOT_AVAILABLE_STRING;
        }
    }
}
