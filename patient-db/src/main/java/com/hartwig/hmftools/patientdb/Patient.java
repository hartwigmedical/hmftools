package com.hartwig.hmftools.patientdb;

import java.util.List;
import java.util.Optional;

import org.jetbrains.annotations.NotNull;

class Patient {
    private final PatientInfo patientInfo;
    private final Optional<TumorData> tumorData;
    private final Optional<List<SystemicTherapyData>> systemicTherapies;
    private final Optional<List<RadioTherapyData>> radioTherapies;
    private final Optional<TreatmentData> treatmentData;

    Patient(@NotNull PatientInfo patientInfo, @NotNull Optional<TumorData> tumorData,
            @NotNull Optional<List<SystemicTherapyData>> systemicTherapies,
            @NotNull Optional<List<RadioTherapyData>> radioTherapies, @NotNull Optional<TreatmentData> treatmentData) {
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
