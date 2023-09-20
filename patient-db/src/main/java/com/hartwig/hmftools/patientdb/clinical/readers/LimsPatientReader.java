package com.hartwig.hmftools.patientdb.clinical.readers;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.patientdb.clinical.curators.DoidNodesResolver;
import com.hartwig.hmftools.patientdb.clinical.datamodel.*;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;
import org.jetbrains.annotations.NotNull;

import java.util.List;

public class LimsPatientReader {

    private static final String DOID_STRING_DELIMITER = ",";

    @NotNull
    private final DoidNodesResolver doidNodesResolver;

    public LimsPatientReader(@NotNull DoidNodesResolver doidNodesResolver) {
        this.doidNodesResolver = doidNodesResolver;
    }

    @NotNull
    public Patient read(@NotNull String patientIdentifier, @NotNull SampleData chosenSample, @NotNull List<SampleData> sequencedSamples) {
        List<DoidNode> resolvedDoidNodes = Lists.newArrayList();
        if (chosenSample.limsTumorDoids() != null){
            resolvedDoidNodes = doidNodesResolver.resolveDoidNodes(List.of(chosenSample.limsTumorDoids().split(DOID_STRING_DELIMITER)));
        }

        CuratedPrimaryTumor curatedPrimaryTumor = ImmutableCuratedPrimaryTumor.builder()
                .location(chosenSample.limsTumorLocation())
                .type(chosenSample.limsTumorType())
                .extraDetails(chosenSample.limsTumorExtra())
                .doidNodes(resolvedDoidNodes)
                .isOverridden(false)
                .build();

        return new Patient(patientIdentifier,
                toBaselineData(curatedPrimaryTumor),
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
    private static BaselineData toBaselineData(@NotNull CuratedPrimaryTumor curatedPrimaryTumor) {
        return ImmutableBaselineData.builder()
                .registrationDate(null)
                .informedConsentDate(null)
                .pifVersion(null)
                .inDatabase(null)
                .outsideEU(null)
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
