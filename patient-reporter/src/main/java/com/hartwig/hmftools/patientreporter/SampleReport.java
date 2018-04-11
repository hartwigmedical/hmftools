package com.hartwig.hmftools.patientreporter;

import java.time.LocalDate;

import com.hartwig.hmftools.patientreporter.util.PatientReportFormat;

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
    public abstract String primaryTumorLocation();

    @Nullable
    public abstract Double pathologyTumorPercentage();

    public String tumorPercentageString() {
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
