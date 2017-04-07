package com.hartwig.hmftools.patientdb.data;

import java.util.List;
import java.util.Optional;

import org.jetbrains.annotations.NotNull;

public class Patient {
    @NotNull
    private final PatientInfo patientInfo;
    @NotNull
    private final Optional<TumorData> tumorData;
    @NotNull
    private final List<SystemicTherapyData> systemicTherapies;
    @NotNull
    private final List<RadioTherapyData> radioTherapies;
    @NotNull
    private final Optional<TreatmentData> treatmentData;
    @NotNull
    private final List<SomaticVariantData> somaticVariants;

    public Patient(@NotNull final PatientInfo patientInfo, @NotNull final Optional<TumorData> tumorData,
            @NotNull final List<SystemicTherapyData> systemicTherapies,
            @NotNull final List<RadioTherapyData> radioTherapies, @NotNull final Optional<TreatmentData> treatmentData,
            @NotNull final List<SomaticVariantData> somaticVariants) {
        this.patientInfo = patientInfo;
        this.tumorData = tumorData;
        this.systemicTherapies = systemicTherapies;
        this.radioTherapies = radioTherapies;
        this.treatmentData = treatmentData;
        this.somaticVariants = somaticVariants;
    }

    @NotNull
    public List<RadioTherapyData> radioTherapies() {
        return radioTherapies;
    }

    @NotNull
    public List<SystemicTherapyData> systemicTherapies() {
        return systemicTherapies;
    }

    @NotNull
    public Optional<TreatmentData> treatmentData() {
        return treatmentData;
    }

    @NotNull
    public Optional<TumorData> tumorData() {
        return tumorData;
    }

    @NotNull
    public PatientInfo patientInfo() {
        return patientInfo;
    }

    @NotNull
    public List<SomaticVariantData> somaticVariants() {
        return somaticVariants;
    }
}
