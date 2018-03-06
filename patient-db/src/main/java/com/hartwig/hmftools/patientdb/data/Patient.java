package com.hartwig.hmftools.patientdb.data;

import java.util.List;

import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;

import org.jetbrains.annotations.NotNull;

public class Patient {
    @NotNull
    private final PatientData patientData;
    @NotNull
    private final PreTreatmentData preTreatmentData;
    @NotNull
    private final List<SampleData> sequencedBiopsies;
    @NotNull
    private final List<BiopsyData> clinicalBiopsies;
    @NotNull
    private final List<BiopsyTreatmentData> treatments;
    @NotNull
    private final List<BiopsyTreatmentResponseData> treatmentResponses;
    @NotNull
    private final List<ValidationFinding> matchFindings;

    public Patient(@NotNull final PatientData patientData, @NotNull final PreTreatmentData preTreatmentData,
            @NotNull final List<SampleData> sequencedBiopsies, @NotNull final List<BiopsyData> clinicalBiopsies,
            @NotNull final List<BiopsyTreatmentData> treatments, @NotNull final List<BiopsyTreatmentResponseData> treatmentResponses,
            @NotNull final List<ValidationFinding> matchFindings) {
        this.patientData = patientData;
        this.preTreatmentData = preTreatmentData;
        this.sequencedBiopsies = sequencedBiopsies;
        this.clinicalBiopsies = clinicalBiopsies;
        this.treatments = treatments;
        this.treatmentResponses = treatmentResponses;
        this.matchFindings = matchFindings;
    }

    @NotNull
    public PatientData patientData() {
        return patientData;
    }

    @NotNull
    public PreTreatmentData preTreatmentData() {
        return preTreatmentData;
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

    @NotNull
    public List<ValidationFinding> matchFindings() {
        return matchFindings;
    }
}
