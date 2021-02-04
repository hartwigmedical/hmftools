package com.hartwig.hmftools.patientdb.clinical.readers.wide;

import java.time.LocalDate;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class WideFiveDays implements WideClinicalData {

    @NotNull
    @Override
    public abstract String widePatientId();

    public abstract boolean dataIsAvailable();

    // When data is not available, all other fields will be null
    @Nullable
    public abstract LocalDate informedConsentDate();

    @Nullable
    public abstract String gender();

    @Nullable
    public abstract Integer birthYear();

    @Nullable
    public abstract LocalDate biopsyDate();

    @Nullable
    public abstract String biopsySite();

    @Nullable
    public abstract String sampleTissue();

    @Nullable
    public abstract String sampleType();

    @Nullable
    public abstract String studyCode();

    @Nullable
    public abstract Boolean participatesInOtherTrials();

    @Nullable
    public abstract String otherTrialCodes();

    @Nullable
    public abstract String otherTrialStartDates();
}
