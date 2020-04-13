package com.hartwig.hmftools.patientdb.readers.wide;

import java.time.LocalDate;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class WideFiveDays {

    @NotNull
    public abstract String patientId();

    public abstract boolean dataIsAvailable();

    @Nullable
    public abstract LocalDate informedConsentDate();

    @Nullable
    public abstract String gender();

    @Nullable
    public abstract Integer birthYear();

    @Nullable
    public abstract LocalDate biopsyDate();

    @NotNull
    public abstract String biopsySite();

    @NotNull
    public abstract String sampleTissue();

    @NotNull
    public abstract String sampleType();

    @NotNull
    public abstract String studyCode();

    @Nullable
    public abstract Boolean participatesInOtherTrials();

    @NotNull
    public abstract String otherTrialCodes();

    @NotNull
    public abstract String otherTrialStartDates();
}
