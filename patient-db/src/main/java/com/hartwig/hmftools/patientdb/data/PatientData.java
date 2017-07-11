package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PatientData {
    @Nullable
    public abstract String cpctId();

    @Nullable
    public abstract LocalDate registrationDate();

    @Nullable
    public abstract String gender();

    @Nullable
    public abstract String ethnicity();

    @Nullable
    public abstract String hospital();

    @Nullable
    public abstract Integer birthYear();

    @Nullable
    public abstract String primaryTumorLocation();

    @Nullable
    public abstract LocalDate deathDate();

    @Nullable
    public abstract String demographyStatus();

    @Nullable
    public abstract String demographyLocked();

    @Nullable
    public abstract String primaryTumorStatus();

    @Nullable
    public abstract String primaryTumorLocked();

    @Nullable
    public abstract String eligibilityStatus();

    @Nullable
    public abstract String eligibilityLocked();

    @Nullable
    public abstract String selectionCriteriaStatus();

    @Nullable
    public abstract String selectionCriteriaLocked();

    @Nullable
    public abstract String deathStatus();

    @Nullable
    public abstract String deathLocked();
}
