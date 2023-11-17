package com.hartwig.hmftools.patientdb.clinical.readers.wide.V2;

import java.time.LocalDate;
import java.util.Optional;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class })
public abstract class PatientData
{
    @NotNull
    public abstract String subjectKey();

    public abstract Optional<LocalDate> informedConsentDate();

    public abstract Optional<LocalDate> registrationDate();

    public abstract Optional<Integer> yearOfBirth();

    public abstract Optional<Gender> gender();

    public abstract Optional<Boolean> specimensCanBeStoredIndefinitelyAndBeUsedForFutureResearch();

    public abstract Optional<Boolean> isDataAvailable();

    public abstract Optional<Boolean> otherTrial(); // boolean doesn't really make sense with this naming.

    public abstract Optional<String> otherTrialCode();

    public abstract Optional<LocalDate> otherTrialDate();

    public abstract Optional<Boolean> hadPreviousChemoTherapy();

    public abstract Optional<String> chemoEstimatedDate(); //TODO what? Why is this a (boolean) string according to specs

    public abstract Optional<LocalDate> chemoLastDoseDate();

    public abstract Optional<Boolean> hadPreviousRadioTherapy();

    public abstract Optional<String> radioEstimatedDate(); //TODO what? Why is this a (boolean) string according to specs

    public abstract Optional<LocalDate> radioLastDoseDate();

    public static ImmutablePatientData.Builder builder()
    {
        return ImmutablePatientData.builder();
    }

    public enum Gender
    {
        MALE,
        FEMALE
    }
}
