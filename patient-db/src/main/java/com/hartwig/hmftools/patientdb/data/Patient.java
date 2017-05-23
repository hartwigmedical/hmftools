package com.hartwig.hmftools.patientdb.data;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class Patient {
    @NotNull
    private final PatientData patientData;
    @NotNull
    private final List<SampleData> sequencedBiopsies;
    @NotNull
    private final List<BiopsyData> clinicalBiopsies;
    @NotNull
    private final List<BiopsyTreatmentData> treatments;
    @NotNull
    private final List<BiopsyTreatmentResponseData> treatmentResponses;

    public Patient(@NotNull final PatientData patientData, @NotNull final List<SampleData> sequencedBiopsies,
            @NotNull final List<BiopsyData> clinicalBiopsies, @NotNull final List<BiopsyTreatmentData> treatments,
            @NotNull final List<BiopsyTreatmentResponseData> treatmentResponses) {
        this.patientData = patientData;
        this.sequencedBiopsies = sequencedBiopsies;
        this.clinicalBiopsies = clinicalBiopsies;
        this.treatments = treatments;
        this.treatmentResponses = treatmentResponses;
    }

    @NotNull
    public PatientData patientInfo() {
        return patientData;
    }

    @NotNull
    public List<SampleData> sequencedBiopsies() {
        return sequencedBiopsies;
    }

    @NotNull
    public List<BiopsyData> clinicalBiopsies() {
        return clinicalBiopsies;
    }

    @NotNull
    public List<BiopsyTreatmentData> treatments() {
        return treatments;
    }

    @NotNull
    public List<BiopsyTreatmentResponseData> treatmentResponses() {
        return treatmentResponses;
    }
}
