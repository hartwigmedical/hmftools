package com.hartwig.hmftools.patientdb.readers.wide;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class WidePreTreatmentData {

    @NotNull
    public abstract String patientId();

    @NotNull
    public abstract String previousTherapy();

    @NotNull
    public abstract String drug1();

    @NotNull
    public abstract String drug2();

    @NotNull
    public abstract String drug3();

    @NotNull
    public abstract String drug4();

    @NotNull
    public abstract String dateLastSystemicTherapy();
}