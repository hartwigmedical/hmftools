package com.hartwig.hmftools.patientdb.clinical.readers;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.clinical.consents.ConsentConfig;
import com.hartwig.hmftools.patientdb.clinical.curators.PrimaryTumorCurator;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BaselineData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.CuratedPrimaryTumor;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableBaselineData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutablePreTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.Patient;
import com.hartwig.hmftools.patientdb.clinical.datamodel.PreTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.SampleData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CorePatientReader {

    @NotNull
    private final PrimaryTumorCurator primaryTumorCurator;

    public CorePatientReader(@NotNull final PrimaryTumorCurator primaryTumorCurator) {
        this.primaryTumorCurator = primaryTumorCurator;
    }

    @NotNull
    public Patient read(@NotNull String patientIdentifier, @Nullable String limsPrimaryTumorLocation,
            @NotNull List<SampleData> sequencedSamples, @NotNull Map<String, ConsentConfig> consentConfigMap, @NotNull String cohort) {
        CuratedPrimaryTumor curatedPrimaryTumor = primaryTumorCurator.search(patientIdentifier, limsPrimaryTumorLocation);
        ConsentConfig extractConsentConfigInfo = consentConfigMap.get(cohort);

        return new Patient(patientIdentifier,
                toBaselineData(curatedPrimaryTumor, extractConsentConfigInfo),
                noPreTreatmentData(),
                sequencedSamples,
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList(),
                Lists.newArrayList());
    }

    @NotNull
    private BaselineData toBaselineData(@NotNull CuratedPrimaryTumor curatedPrimaryTumor,
            @Nullable ConsentConfig extractConsentConfigInfo) {

        return ImmutableBaselineData.builder()
                .registrationDate(null)
                .informedConsentDate(null)
                .pifVersion(extractConsentConfigInfo != null ? extractConsentConfigInfo.pifVersion() : null)
                .inDatabase(extractConsentConfigInfo != null ? extractConsentConfigInfo.inHMF() : null)
                .outsideEU(extractConsentConfigInfo != null ? extractConsentConfigInfo.outsideEU() : null)
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
