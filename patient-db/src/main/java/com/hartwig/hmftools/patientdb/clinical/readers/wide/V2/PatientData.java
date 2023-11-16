package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import java.time.LocalDate;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class PatientData
{
    @NotNull
    public abstract String subjectKey();

    public abstract LocalDate informedConsentDate();

    public abstract LocalDate registrationDate();

    public abstract int yearOfBirth();

    public abstract Gender gender();

    public abstract boolean specimensCanBeStoredIndefinitelyAndBeUsedForFutureResearch();

    public abstract boolean dataAvailable();

    public abstract boolean otherTrial(); // boolean doesn't really make sense with this naming.

    public abstract String otherTrialCode();

    public abstract LocalDate otherTrialDate();

    public abstract boolean hadPreviousChemoTherapy();

    public abstract String chemoEstimatedDate(); //TODO what? Why is this a (boolean) string according to specs

    public abstract LocalDate chemoDateLastDose();

    public abstract boolean hadPreviousRadioTherapy();

    public abstract String radioEstimatedDate(); //TODO what? Why is this a (boolean) string according to specs

    public abstract LocalDate radioDateLastDose();

    public enum Gender
    {
        MALE,
        FEMALE
    }
}
