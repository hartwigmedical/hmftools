package com.hartwig.hmftools.patientdb.clinical.datamodel;

import java.time.LocalDate;

import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class DatamodelTestFactory {

    private DatamodelTestFactory() {
    }

    @NotNull
    public static ImmutableSampleData.Builder sampleBuilder(@NotNull LocalDate arrivalDate) {
        return ImmutableSampleData.builder()
                .sampleId("sample-" + arrivalDate.toString())
                .sampleBarcode("ABC")
                .cohortId(Strings.EMPTY)
                .sequenced(false)
                .isSomaticTumorSample(true)
                .requiresCuratedPrimaryTumor(true)
                .setName(Strings.EMPTY)
                .arrivalDate(arrivalDate)
                .pathologyTumorPercentage("N/A")
                .pathologySampleId("N/A");
    }

    @NotNull
    public static ImmutableBaselineData.Builder baselineBuilder() {
        return ImmutableBaselineData.builder()
                .curatedPrimaryTumor(ImmutableCuratedPrimaryTumor.builder().searchTerm(Strings.EMPTY).isOverridden(false).build())
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
