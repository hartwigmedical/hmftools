package com.hartwig.hmftools.patientdb.readers.wide;

import java.time.LocalDate;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class WideResponseData {

    @NotNull
    public abstract String patientId();

    public abstract int timePoint();

    @NotNull
    public abstract LocalDate date();

    public abstract boolean recistDone();

    @NotNull
    public abstract String recistResponse();

    @NotNull
    public abstract String noRecistResponse();

    @NotNull
    public abstract String noRecistReasonStopTreatment();

    @NotNull
    public abstract String noRecistReasonStopTreatmentOther();

}
