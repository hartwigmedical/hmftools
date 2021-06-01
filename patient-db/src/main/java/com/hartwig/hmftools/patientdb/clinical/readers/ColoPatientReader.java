package com.hartwig.hmftools.patientdb.clinical.readers;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BaselineData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.CuratedPrimaryTumor;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableBaselineData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableCuratedPrimaryTumor;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutablePreTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.Patient;
import com.hartwig.hmftools.patientdb.clinical.datamodel.PreTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ColoPatientReader {

    public ColoPatientReader() {
    }

    @NotNull
    public Patient read(@NotNull String coloSampleId) {
        boolean isColo829 = coloSampleId.equals("COLO829T");

        String patientId = isColo829 ? "COLO829" : Strings.EMPTY;

        ImmutableCuratedPrimaryTumor.Builder primaryTumorBuilder =
                ImmutableCuratedPrimaryTumor.builder().searchTerm(Strings.EMPTY).isOverridden(false);

        CuratedPrimaryTumor curatedPrimaryTumor = isColo829 ? primaryTumorBuilder.location("Skin")
                .subLocation(Strings.EMPTY)
                .type("Melanoma")
                .subType(Strings.EMPTY)
                .extraDetails(Strings.EMPTY)
                .build() : primaryTumorBuilder.build();

        return new Patient(patientId,
                toBaselineData(curatedPrimaryTumor),
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
    private static BaselineData toBaselineData(@NotNull CuratedPrimaryTumor curatedPrimaryTumor) {
        return ImmutableBaselineData.builder()
                .registrationDate(null)
                .informedConsentDate(null)
                .gender(null)
                .hospital(null)
                .birthYear(null)
                .curatedPrimaryTumor(curatedPrimaryTumor)
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
