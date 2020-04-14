package com.hartwig.hmftools.patientdb.readers.wide;

import java.time.LocalDate;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class WideBiopsyData {

    @NotNull
    public abstract String patientId();

    @NotNull
    public abstract String pathologySampleId();

    @Nullable
    public abstract LocalDate biopsyDate();

    @Nullable
    public abstract Boolean hasReceivedSuccessfulReport();
}