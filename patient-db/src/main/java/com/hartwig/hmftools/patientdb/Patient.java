package com.hartwig.hmftools.patientdb;

import java.util.List;
import java.util.Optional;

import org.jetbrains.annotations.NotNull;

class Patient {
    @NotNull
    private final PatientInfo patientInfo;
    @NotNull
    private final Optional<TumorData> tumorData;
    @NotNull
    private final Optional<List<SystemicTherapyData>> systemicTherapies;
    @NotNull
    private final Optional<List<RadioTherapyData>> radioTherapies;
    @NotNull
    private final Optional<TreatmentData> treatmentData;

    Patient(@NotNull final PatientInfo patientInfo, @NotNull final Optional<TumorData> tumorData,
            @NotNull final Optional<List<SystemicTherapyData>> systemicTherapies,
            @NotNull final Optional<List<RadioTherapyData>> radioTherapies,
            @NotNull final Optional<TreatmentData> treatmentData) {
        this.patientInfo = patientInfo;
        this.tumorData = tumorData;
        this.systemicTherapies = systemicTherapies;
        this.radioTherapies = radioTherapies;
        this.treatmentData = treatmentData;
    }

    Optional<List<RadioTherapyData>> radioTherapies() {
        return radioTherapies;
    }

    Optional<List<SystemicTherapyData>> systemicTherapies() {
        return systemicTherapies;
    }

    Optional<TreatmentData> treatmentData() {
        return treatmentData;
    }

    Optional<TumorData> tumorData() {
        return tumorData;
    }

    PatientInfo patientInfo() {
        return patientInfo;
    }
}
