package com.hartwig.hmftools.patientreporter;

import java.time.LocalDate;

import com.google.gson.annotations.SerializedName;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.lims.LimsSampleType;

import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class SampleReport {

    @NotNull
    public abstract String sampleId();

    @NotNull
    public abstract String barcodeTumor();

    @NotNull
    public abstract String barcodeReference();

    @Nullable
    public abstract PatientTumorLocation patientTumorLocation();

    @NotNull
    public abstract String purityShallowSeq();

    @NotNull
    public abstract String pathologyTumorPercentage();

    @Nullable
    public abstract LocalDate tumorArrivalDate();

    @Nullable
    public abstract LocalDate bloodArrivalDate();

    @NotNull
    public abstract String labProcedures();

    @Nullable
    public abstract String recipient();

    @NotNull
    public abstract String projectName();

    @NotNull
    public abstract String contactNames();

    @NotNull
    public abstract String contactEmails();

    @NotNull
    public abstract String submission();

    @Nullable
    public abstract String hospitalPatientId();

    @NotNull
    @Value.Derived
    public String primaryTumorLocationString() {
        PatientTumorLocation type = patientTumorLocation();
        return type != null ? type.primaryTumorLocation() : Strings.EMPTY;
    }

    @NotNull
    @Value.Derived
    public String cancerSubTypeString() {
        PatientTumorLocation type = patientTumorLocation();
        return type != null ? type.cancerSubtype() : Strings.EMPTY;
    }

    @NotNull
    @Value.Derived
    public String buildReportTitle(@NotNull String title) {
        LimsSampleType type = LimsSampleType.fromSampleId(sampleId());

        String patientNumber = hospitalPatientId();
        if (type == LimsSampleType.CORE && patientNumber == null) {
            throw new IllegalStateException("CORE sample present without patient number: " + sampleId());
        }

        return type == LimsSampleType.CORE ? title + " - " + sampleId() + " (" + hospitalPatientId() + ")" : title + " - " + sampleId();
    }
}
