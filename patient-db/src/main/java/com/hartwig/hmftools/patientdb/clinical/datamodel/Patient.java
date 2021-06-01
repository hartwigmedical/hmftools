package com.hartwig.hmftools.patientdb.clinical.datamodel;

import java.util.List;

import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.ValidationFinding;

import org.jetbrains.annotations.NotNull;

public class Patient {

    @NotNull
    private final String patientIdentifier;
    @NotNull
    private final BaselineData baselineData;
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
    private final List<TumorMarkerData> tumorMarkers;
    @NotNull
    private final List<RanoMeasurementData> ranoMeasurements;
    @NotNull
    private final List<ValidationFinding> matchFindings;

    public Patient(@NotNull final String patientIdentifier, @NotNull final BaselineData baselineData,
            @NotNull final PreTreatmentData preTreatmentData, @NotNull final List<SampleData> sequencedBiopsies,
            @NotNull final List<BiopsyData> clinicalBiopsies, @NotNull final List<BiopsyTreatmentData> treatments,
            @NotNull final List<BiopsyTreatmentResponseData> treatmentResponses, @NotNull final List<TumorMarkerData> tumorMarkers,
            @NotNull final List<RanoMeasurementData> ranoMeasurements, @NotNull final List<ValidationFinding> matchFindings) {
        this.patientIdentifier = patientIdentifier;
        this.baselineData = baselineData;
        this.preTreatmentData = preTreatmentData;
        this.sequencedBiopsies = sequencedBiopsies;
        this.clinicalBiopsies = clinicalBiopsies;
        this.treatments = treatments;
        this.treatmentResponses = treatmentResponses;
        this.tumorMarkers = tumorMarkers;
        this.ranoMeasurements = ranoMeasurements;
        this.matchFindings = matchFindings;
    }

    @NotNull
    public String patientIdentifier() {
        return patientIdentifier;
    }

    @NotNull
    public BaselineData baselineData() {
        return baselineData;
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
    public List<TumorMarkerData> tumorMarkers() {
        return tumorMarkers;
    }

    @NotNull
    public List<RanoMeasurementData> ranoMeasurements() {
        return ranoMeasurements;
    }

    @NotNull
    public List<ValidationFinding> matchFindings() {
        return matchFindings;
    }
}
