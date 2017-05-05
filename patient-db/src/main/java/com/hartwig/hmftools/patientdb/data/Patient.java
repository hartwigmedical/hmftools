package com.hartwig.hmftools.patientdb.data;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class Patient {
    @NotNull
    private final PatientInfo patientInfo;
    @NotNull
    private final List<BiopsyLimsData> sequencedBiopsies;
    @NotNull
    private final List<BiopsyTreatmentData> treatments;
    @NotNull
    private final List<BiopsyTreatmentResponseData> treatmentResponses;
    @NotNull
    private final List<BiopsyClinicalData> clinicalBiopsies;
    @NotNull
    private final List<SomaticVariantData> somaticVariants;

    public Patient(@NotNull final PatientInfo patientInfo, @NotNull final List<BiopsyLimsData> sequencedBiopsies,
            @NotNull final List<BiopsyTreatmentData> treatments,
            @NotNull final List<BiopsyTreatmentResponseData> treatmentResponses,
            @NotNull final List<BiopsyClinicalData> clinicalBiopsies,
            @NotNull final List<SomaticVariantData> somaticVariants) {
        this.patientInfo = patientInfo;
        this.sequencedBiopsies = sequencedBiopsies;
        this.treatments = treatments;
        this.treatmentResponses = treatmentResponses;
        this.clinicalBiopsies = clinicalBiopsies;
        this.somaticVariants = somaticVariants;
    }

    @NotNull
    public PatientInfo patientInfo() {
        return patientInfo;
    }

    @NotNull
    public List<SomaticVariantData> somaticVariants() {
        return somaticVariants;
    }

    @NotNull
    public List<BiopsyLimsData> sequencedBiopsies() {
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
    public List<BiopsyClinicalData> clinicalBiopsies() {
        return clinicalBiopsies;
    }

}
