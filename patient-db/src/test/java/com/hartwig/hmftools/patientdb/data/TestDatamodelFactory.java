package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;

import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;

import org.jetbrains.annotations.NotNull;

public final class TestDatamodelFactory {

    private TestDatamodelFactory() {
    }

    @NotNull
    public static ImmutableSampleData.Builder sampleBuilder(@NotNull LocalDate arrivalDate) {
        return ImmutableSampleData.builder().sampleId("sample-" + arrivalDate.toString()).arrivalDate(arrivalDate);
    }

    @NotNull
    public static ImmutableBaselineData.Builder baseLineBuilder() {
        return ImmutableBaselineData.builder()
                .cancerType(ImmutableCuratedCancerType.of(null, null, null))
                .demographyLocked(false)
                .demographyStatus(FormStatusState.UNKNOWN)
                .primaryTumorLocked(false)
                .primaryTumorStatus(FormStatusState.UNKNOWN)
                .selectionCriteriaLocked(false)
                .selectionCriteriaStatus(FormStatusState.UNKNOWN)
                .eligibilityLocked(false)
                .eligibilityStatus(FormStatusState.UNKNOWN)
                .informedConsentLocked(false)
                .informedConsentStatus(FormStatusState.UNKNOWN)
                .deathLocked(false)
                .deathStatus(FormStatusState.UNKNOWN);
    }

    @NotNull
    public static ImmutableBiopsyData.Builder biopsyBuilder() {
        return ImmutableBiopsyData.builder().id(1).formStatus(FormStatusState.UNKNOWN).formLocked(false);
    }

    @NotNull
    public static ImmutableBiopsyTreatmentData.Builder biopsyTreatmentBuilder() {
        return ImmutableBiopsyTreatmentData.builder().id(1).formStatus(FormStatusState.UNKNOWN).formLocked(false);
    }

    @NotNull
    public static ImmutableDrugData.Builder drugBuilder() {
        return ImmutableDrugData.builder();
    }

    @NotNull
    public static ImmutableBiopsyTreatmentResponseData.Builder biopsyTreatmentResponseBuilder() {
        return ImmutableBiopsyTreatmentResponseData.builder().formStatus(FormStatusState.UNKNOWN).formLocked(false);
    }
}
