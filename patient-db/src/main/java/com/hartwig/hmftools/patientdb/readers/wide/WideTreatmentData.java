package com.hartwig.hmftools.patientdb.readers.wide;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class WideTreatmentData {

    @NotNull
    public abstract String sampleId();

    @NotNull
    public abstract String drugCode();

    @NotNull
    public abstract String drug();

    @NotNull
    public abstract String startDate();

    @NotNull
    public abstract String endDate();
}