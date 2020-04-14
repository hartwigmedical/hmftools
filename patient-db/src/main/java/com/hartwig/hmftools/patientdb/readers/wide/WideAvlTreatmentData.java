package com.hartwig.hmftools.patientdb.readers.wide;

import java.time.LocalDate;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class WideAvlTreatmentData {

    @NotNull
    public abstract String patientId();

    @NotNull
    public abstract String drugCode();

    @NotNull
    public abstract String drug();

    @Nullable
    public abstract LocalDate startDate();

    @Nullable
    public abstract LocalDate endDate();
}