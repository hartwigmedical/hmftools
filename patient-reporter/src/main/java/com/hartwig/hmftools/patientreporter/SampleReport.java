package com.hartwig.hmftools.patientreporter;

import java.time.LocalDate;

import com.hartwig.hmftools.common.ecrf.projections.PatientCancerType;
import com.hartwig.hmftools.patientreporter.util.PatientReportFormat;

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
    public abstract PatientCancerType patientCancerType();

    @Nullable
    public abstract Double pathologyTumorPercentage();

    @NotNull
    @Value.Derived
    public String primaryTumorLocationString() {
        PatientCancerType type = patientCancerType();
        return type != null ? type.primaryTumorLocation() : Strings.EMPTY;
    }

    @NotNull
    @Value.Derived
    public String cancerSubTypeString() {
        PatientCancerType type = patientCancerType();
        return type != null ? type.cancerSubtype() : Strings.EMPTY;
    }

    @NotNull
    @Value.Derived
    public String pathologyTumorPercentageString() {
        return PatientReportFormat.formatNullablePercent(pathologyTumorPercentage());
    }

    @Nullable
    public abstract LocalDate tumorArrivalDate();

    @Nullable
    public abstract LocalDate bloodArrivalDate();

    @NotNull
    public abstract String labProcedures();

    @Nullable
    public abstract String recipient();
}
