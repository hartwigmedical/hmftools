package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;

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
    public abstract String gender();

    @Nullable
    public abstract String hospital();

    @Nullable
    public abstract Integer birthYear();

    @NotNull
    @Value.Default
    public CuratedTumorLocation primaryTumorLocation() {
        return ImmutableCuratedTumorLocation.of(null, null, null);
    }

    @Nullable
    public abstract LocalDate deathDate();

    @NotNull
    @Value.Default
    public String demographyStatus() {
        return "";
    }

    @NotNull
    @Value.Default
    public String demographyLocked() {
        return "";
    }

    @NotNull
    @Value.Default
    public String primaryTumorStatus() {
        return "";
    }

    @NotNull
    @Value.Default
    public String primaryTumorLocked() {
        return "";
    }

    @NotNull
    @Value.Default
    public String eligibilityStatus() {
        return "";
    }

    @NotNull
    @Value.Default
    public String eligibilityLocked() {
        return "";
    }

    @NotNull
    @Value.Default
    public String selectionCriteriaStatus() {
        return "";
    }

    @NotNull
    @Value.Default
    public String selectionCriteriaLocked() {
        return "";
    }

    @NotNull
    @Value.Default
    public String deathStatus() {
        return "";
    }

    @NotNull
    @Value.Default
    public String deathLocked() {
        return "";
    }
}
