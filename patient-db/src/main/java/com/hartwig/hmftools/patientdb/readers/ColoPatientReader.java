package com.hartwig.hmftools.patientdb.readers;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.data.BaselineData;
import com.hartwig.hmftools.patientdb.data.CuratedTumorLocationV2;
import com.hartwig.hmftools.patientdb.data.ImmutableBaselineData;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedTumorLocationV2;
import com.hartwig.hmftools.patientdb.data.ImmutablePreTreatmentData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PreTreatmentData;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ColoPatientReader {

    @NotNull
    public Patient read(@NotNull String coloSampleId) {
        boolean isColo829 = coloSampleId.equals("COLO829T");

        String patientId = isColo829 ? "COLO829" : Strings.EMPTY;

        CuratedTumorLocationV2 curatedTumorLocationV2 = isColo829
                ? ImmutableCuratedTumorLocationV2.builder().primaryTumorLocation("Skin").primaryTumorType("Melanoma").build()
                : ImmutableCuratedTumorLocationV2.builder().build();

        return new Patient(patientId,
                toBaselineData(curatedTumorLocationV2),
                noPreTreatmentData(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList());
    }

    @NotNull
    private static BaselineData toBaselineData(
            @NotNull CuratedTumorLocationV2 curatedTumorLocationV2) {
        return ImmutableBaselineData.builder()
                .registrationDate(null)
                .informedConsentDate(null)
                .gender(null)
                .hospital(null)
                .birthYear(null)
                .curatedTumorLocationV2(curatedTumorLocationV2)
                .deathDate(null)
                .demographyStatus(FormStatus.undefined())
                .primaryTumorStatus(FormStatus.undefined())
                .informedConsentStatus(FormStatus.undefined())
                .eligibilityStatus(FormStatus.undefined())
                .selectionCriteriaStatus(FormStatus.undefined())
                .deathStatus(FormStatus.undefined())
                .build();
    }

    @NotNull
    private static PreTreatmentData noPreTreatmentData() {
        return ImmutablePreTreatmentData.builder()
                .treatmentGiven(null)
                .radiotherapyGiven(null)
                .drugs(Lists.newArrayList())
                .formStatus(FormStatus.undefined())
                .build();
    }
}
