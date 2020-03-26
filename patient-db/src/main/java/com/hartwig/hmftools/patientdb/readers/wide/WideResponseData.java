package com.hartwig.hmftools.patientdb.readers.wide;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class WideResponseData {

    @NotNull
    public abstract String patientId();

    @NotNull
    public abstract String timePoint();

    @NotNull
    public abstract String date();

    @NotNull
    public abstract String recistNotDone();

    @NotNull
    public abstract String responseAccordingRecist();

    @NotNull
    public abstract String clinicalDecision();

    @NotNull
    public abstract String reasonStopTreatment();

    @NotNull
    public abstract String reasonStopTreatmentOther();

}
