package com.hartwig.hmftools.patientdb.clinical.datamodel;

import java.time.LocalDate;

import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class RanoMeasurementData {

    @NotNull
    public abstract String patientId();

    @Nullable
    public abstract String therapyGiven();

    @Nullable
    public abstract LocalDate responseDate();

    @Nullable
    public abstract String targetLesionResponse();

    @Nullable
    public abstract String noTargetLesionResponse();

    @Nullable
    public abstract String overallResponse();

    @NotNull
    public abstract FormStatus formStatus();
}
