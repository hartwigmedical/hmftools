package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;

import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PatientData {
    @NotNull
    public abstract String cpctId();

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
    @Value.Default
    public CuratedCancerType cancerType() {
        return ImmutableCuratedCancerType.of(null, null, null);
    }

    @Nullable
    public abstract LocalDate deathDate();

    @NotNull
    @Value.Default
    public FormStatusState demographyStatus() {
        return FormStatusState.UNKNOWN;
    }

    @Value.Default
    public boolean demographyLocked() {
        return false;
    }

    @NotNull
    @Value.Default
    public FormStatusState primaryTumorStatus() {
        return FormStatusState.UNKNOWN;
    }

    @Value.Default
    public boolean primaryTumorLocked() {
        return false;
    }

    @NotNull
    @Value.Default
    public FormStatusState informedConsentStatus() {
        return FormStatusState.UNKNOWN;
    }

    @Value.Default
    public boolean informedConsentLocked() {
        return false;
    }

    @NotNull
    @Value.Default
    public FormStatusState eligibilityStatus() {
        return FormStatusState.UNKNOWN;
    }

    @Value.Default
    public boolean eligibilityLocked() {
        return false;
    }

    @NotNull
    @Value.Default
    public FormStatusState selectionCriteriaStatus() {
        return FormStatusState.UNKNOWN;
    }

    @Value.Default
    public boolean selectionCriteriaLocked() {
        return false;
    }

    @NotNull
    @Value.Default
    public FormStatusState deathStatus() {
        return FormStatusState.UNKNOWN;
    }

    @Value.Default
    public boolean deathLocked() {
        return false;
    }
}
