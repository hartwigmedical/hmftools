package com.hartwig.hmftools.patientdb.readers.wide;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class WideFiveDays {

    @NotNull
    public abstract String patientId();

    @NotNull
    public abstract String wideId();

    @NotNull
    public abstract String birthYear();

    @NotNull
    public abstract String gender();

    @NotNull
    public abstract String informedConsent();

    @NotNull
    public abstract String usedForFutureResearch();

    @NotNull
    public abstract String dateShared();

    @NotNull
    public abstract String dateBiopsy();

    @NotNull
    public abstract String biopsySite();

    @NotNull
    public abstract String sampleTissue();

    @NotNull
    public abstract String sampleType();

    @NotNull
    public abstract String studyCode();

    @NotNull
    public abstract String otherTrial();

    @NotNull
    public abstract String otherTrialCode();

    @NotNull
    public abstract String dateStartOtherTrial();
}
