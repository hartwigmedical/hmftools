package com.hartwig.hmftools.patientdb.data;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class Patient {
    @NotNull
    private final PatientData patientData;
    @NotNull
    private final List<SampleData> sequencedBiopsies;
    @NotNull
    private final List<BiopsyTreatmentData> treatments;
    @NotNull
    private final List<BiopsyTreatmentResponseData> treatmentResponses;
    @NotNull
    private final List<BiopsyData> clinicalBiopsies;

    public Patient(@NotNull final PatientData patientData, @NotNull final List<SampleData> sequencedBiopsies,
            @NotNull final List<BiopsyTreatmentData> treatments,
            @NotNull final List<BiopsyTreatmentResponseData> treatmentResponses,
            @NotNull final List<BiopsyData> clinicalBiopsies) {
        this.patientData = patientData;
        this.sequencedBiopsies = sequencedBiopsies;
        this.treatments = treatments;
        this.treatmentResponses = treatmentResponses;
        this.clinicalBiopsies = clinicalBiopsies;
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
    public List<BiopsyTreatmentData> treatments() {
        return treatments;
    }

    @NotNull
    public List<BiopsyTreatmentResponseData> treatmentResponses() {
        return treatmentResponses;
    }

    @NotNull
    public List<BiopsyData> clinicalBiopsies() {
        return clinicalBiopsies;
    }

}
