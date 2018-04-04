package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;

import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class BaselineData {
    @Nullable
    public abstract LocalDate registrationDate();

    @Nullable
    public abstract LocalDate informedConsentDate();

    @Nullable
    public abstract String gender();

    @Nullable
    public abstract String hospital();

    @Nullable
    public abstract Integer birthYear();

    @NotNull
    public abstract CuratedCancerType cancerType();

    @Nullable
    public abstract LocalDate deathDate();

    @NotNull
    public abstract FormStatusState demographyStatus();

    public abstract boolean demographyLocked();

    @NotNull
    public abstract FormStatusState primaryTumorStatus();

    public abstract boolean primaryTumorLocked();

    @NotNull
    public abstract FormStatusState informedConsentStatus();

    public abstract boolean informedConsentLocked();

    @NotNull
    public abstract FormStatusState eligibilityStatus();

    public abstract boolean eligibilityLocked();

    @NotNull
    public abstract FormStatusState selectionCriteriaStatus();

    public abstract boolean selectionCriteriaLocked();

    @NotNull
    public abstract FormStatusState deathStatus();

    public abstract boolean deathLocked();
}
