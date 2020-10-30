package com.hartwig.hmftools.patientdb.readers;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.data.BaselineData;
import com.hartwig.hmftools.patientdb.data.CuratedTumorLocation;
import com.hartwig.hmftools.patientdb.data.ImmutableBaselineData;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedTumorLocation;
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

        CuratedTumorLocation curatedTumorLocation = isColo829 ? ImmutableCuratedTumorLocation.builder()
                .primaryTumorLocation("Skin")
                .primaryTumorSubLocation(Strings.EMPTY)
                .primaryTumorType("Melanoma")
                .primaryTumorSubType(Strings.EMPTY)
                .primaryTumorExtraDetails(Strings.EMPTY)
                .build() : ImmutableCuratedTumorLocation.builder().build();

        return new Patient(patientId,
                toBaselineData(curatedTumorLocation),
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
    private static BaselineData toBaselineData(@NotNull CuratedTumorLocation curatedTumorLocation) {
        return ImmutableBaselineData.builder()
                .registrationDate(null)
                .informedConsentDate(null)
                .gender(null)
                .hospital(null)
                .birthYear(null)
                .curatedTumorLocation(curatedTumorLocation)
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
