package com.hartwig.hmftools.patientdb.clinical.datamodel;

import java.time.LocalDate;

import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class TumorMarkerData {

    @NotNull
    public abstract String patientId();

    @Nullable
    public abstract LocalDate date();

    @Nullable
    public abstract String marker();

    @Nullable
    public abstract String measurement();

    @Nullable
    public abstract String unit();

    @NotNull
    public abstract FormStatus formStatus();
}
