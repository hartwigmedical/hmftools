package com.hartwig.hmftools.patientdb.data;

import java.time.LocalDate;

import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;

import org.jetbrains.annotations.NotNull;

public final class TestDatamodelFactory {

    private TestDatamodelFactory() {
    }

    @NotNull
    public static ImmutableSampleData.Builder sampleBuilder(@NotNull LocalDate arrivalDate) {
        return ImmutableSampleData.builder().sampleId("sample-" + arrivalDate.toString()).arrivalDate(arrivalDate);
    }

    @NotNull
    public static ImmutableBaselineData.Builder baselineBuilder() {
        return ImmutableBaselineData.builder().curatedTumorLocation(ImmutableCuratedTumorLocation.of(null, null, null))
                .demographyStatus(FormStatus.undefined())
                .primaryTumorStatus(FormStatus.undefined())
                .selectionCriteriaStatus(FormStatus.undefined())
                .eligibilityStatus(FormStatus.undefined())
                .informedConsentStatus(FormStatus.undefined())
                .deathStatus(FormStatus.undefined());
    }

    @NotNull
    public static ImmutableBiopsyData.Builder biopsyBuilder() {
        return ImmutableBiopsyData.builder().id(1).formStatus(FormStatus.undefined());
    }

    @NotNull
    public static ImmutableBiopsyTreatmentData.Builder biopsyTreatmentBuilder() {
        return ImmutableBiopsyTreatmentData.builder().id(1).formStatus(FormStatus.undefined());
    }

    @NotNull
    public static ImmutableDrugData.Builder drugBuilder() {
        return ImmutableDrugData.builder();
    }

    @NotNull
    public static ImmutableBiopsyTreatmentResponseData.Builder biopsyTreatmentResponseBuilder() {
        return ImmutableBiopsyTreatmentResponseData.builder().formStatus(FormStatus.undefined());
    }
}
