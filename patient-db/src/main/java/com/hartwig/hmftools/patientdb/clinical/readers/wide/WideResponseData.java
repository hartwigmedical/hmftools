package com.hartwig.hmftools.patientdb.clinical.readers.wide;

import java.time.LocalDate;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class WideResponseData implements WideClinicalData {

    @NotNull
    @Override
    public abstract String widePatientId();

    public abstract int timePoint();

    // Response date can be null if no response was taken because patient died before first response
    @Nullable
    public abstract LocalDate date();

    public abstract boolean recistDone();

    @NotNull
    public abstract String recistResponse();

    // When recist is done, these fields will be null.
    @Nullable
    public abstract String noRecistResponse();

    @Nullable
    public abstract String noRecistReasonStopTreatment();

    @Nullable
    public abstract String noRecistReasonStopTreatmentOther();

}
