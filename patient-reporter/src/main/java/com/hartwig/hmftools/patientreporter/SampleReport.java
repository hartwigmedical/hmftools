package com.hartwig.hmftools.patientreporter;

import java.time.LocalDate;

import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;

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

    @Nullable
    public abstract PatientTumorLocation patientTumorLocation();

    @NotNull
    public abstract String refBarcode();

    @Nullable
    public abstract LocalDate refArrivalDate();

    @NotNull
    public abstract String tumorBarcode();

    @Nullable
    public abstract LocalDate tumorArrivalDate();

    @NotNull
    public abstract String purityShallowSeq();

    @NotNull
    public abstract String pathologyTumorPercentage();

    @NotNull
    public abstract String labProcedures();

    @NotNull
    public abstract String projectName();

    @NotNull
    public abstract String requesterName();

    @NotNull
    public abstract String requesterEmail();

    @Nullable
    public abstract String addressee();

    @NotNull
    public abstract String hospitalName();

    @NotNull
    public abstract String hospitalPIName();

    @NotNull
    public abstract String hospitalPIEmail();

    @NotNull
    public abstract String submissionId();

    @NotNull
    public abstract String hospitalPatientId();

    @NotNull
    public abstract String hospitalPathologySampleId();

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
}
