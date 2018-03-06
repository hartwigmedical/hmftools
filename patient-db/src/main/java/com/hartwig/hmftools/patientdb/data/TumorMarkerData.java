package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;

import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;

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
    public abstract FormStatusState formStatus();

    public abstract boolean formLocked();
}
