package com.hartwig.hmftools.patientdb.clinical.datamodel;

import java.time.LocalDate;

import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;

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
    public abstract String pifVersion();

    @Nullable
    public abstract Boolean inDatabase();

    @Nullable
    public abstract Boolean outsideEU();

    @Nullable
    public abstract String gender();

    @Nullable
    public abstract String hospital();

    @Nullable
    public abstract Integer birthYear();

    @NotNull
    public abstract CuratedPrimaryTumor curatedPrimaryTumor();

    @Nullable
    public abstract LocalDate deathDate();

    @NotNull
    public abstract FormStatus demographyStatus();

    @NotNull
    public abstract FormStatus primaryTumorStatus();

    @NotNull
    public abstract FormStatus informedConsentStatus();

    @NotNull
    public abstract FormStatus eligibilityStatus();

    @NotNull
    public abstract FormStatus selectionCriteriaStatus();

    @NotNull
    public abstract FormStatus deathStatus();
}
