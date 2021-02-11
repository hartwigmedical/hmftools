package com.hartwig.hmftools.ckb.interpretation.clinicaltrial;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodelinterpretation.ImmutableCkbEntryInterpretation;
import com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ClinicalTrialLocation;
import com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ImmutableClinicalTrial;
import com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ImmutableClinicalTrialContact;
import com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ImmutableClinicalTrialLocation;
import com.hartwig.hmftools.ckb.datamodelinterpretation.indication.ImmutableIndication;
import com.hartwig.hmftools.ckb.interpretation.common.therapyinterpretation.TherapyInterpretationFactory;
import com.hartwig.hmftools.ckb.interpretation.common.variantinterpretation.VariantInterpretationFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrialContact;
import com.hartwig.hmftools.ckb.json.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.IndicationInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
import com.hartwig.hmftools.ckb.json.indication.Indication;
import com.hartwig.hmftools.ckb.json.molecularprofile.MolecularProfile;
import com.hartwig.hmftools.ckb.json.therapy.Therapy;

import org.jetbrains.annotations.NotNull;

public class ClinicalTrialFactory {

    private ClinicalTrialFactory() {

    }

    public static void interpretClinicalTrials(@NotNull MolecularProfile molecularProfile, @NotNull CkbJsonDatabase ckbEntry,
            @NotNull ImmutableCkbEntryInterpretation.Builder outputBuilder) {
        for (ClinicalTrialInfo clinicalTrialInfo : molecularProfile.variantAssociatedClinicalTrials()) {
            ImmutableClinicalTrialInterpretation.Builder outputBuilderClinicalInterpretation =
                    ImmutableClinicalTrialInterpretation.builder();

            for (ClinicalTrial clinicalTrial : ckbEntry.clinicalTrials()) {
                if (clinicalTrialInfo.nctId().equals(clinicalTrial.nctId())) {
                    outputBuilderClinicalInterpretation.clinicalTrials(ImmutableClinicalTrial.builder()
                            .nctId(clinicalTrial.nctId())
                            .title(clinicalTrial.title())
                            .phase(clinicalTrial.phase())
                            .recruitment(clinicalTrial.recruitment())
                            .ageGroups(clinicalTrial.ageGroups())
                            .gender(clinicalTrial.gender())
                            .variantRequirement(clinicalTrial.variantRequirements())
                            .sponsor(clinicalTrial.sponsors())
                            .updateDate(clinicalTrial.updateDate())
                            .clinicalTrialVariantRequirementDetails(VariantInterpretationFactory.extractProfileName(clinicalTrial.variantRequirementDetails(),
                                    molecularProfile,
                                    ckbEntry))
                            .locations(extractClinicalTrialLocation(clinicalTrial.clinicalTrialLocations()))
                            .build());

                    for (IndicationInfo indicationInfo : clinicalTrial.indications()) {
                        for (Indication indication : ckbEntry.indications()) {
                            if (indicationInfo.id().equals(indication.id())) {
                                outputBuilderClinicalInterpretation.addIndications(ImmutableIndication.builder()
                                        .id(indication.id())
                                        .name(indication.name())
                                        .source(indication.source())
                                        .definition(indication.definition())
                                        .currentPreferredTerm(indication.currentPreferredTerm())
                                        .lastUpdateDateFromDO(indication.lastUpdateDateFromDO())
                                        .altIds(indication.altIds())
                                        .termId(indication.termId())
                                        .build());
                            }
                        }
                    }

                    for (TherapyInfo therapyInfo : clinicalTrial.therapies()) {
                        for (Therapy therapy : ckbEntry.therapies()) {
                            if (therapyInfo.id() == therapy.id()) {

                                outputBuilderClinicalInterpretation.addTherapyInterpretations(TherapyInterpretationFactory.extractTherapyInterpretation(
                                        therapy,
                                        ckbEntry,
                                        molecularProfile));
                            }
                        }
                    }
                }
            }
            outputBuilder.addClinicalTrialInterpretation(outputBuilderClinicalInterpretation.build());
        }
    }

    @NotNull
    private static List<ClinicalTrialLocation> extractClinicalTrialLocation(
            @NotNull List<com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrialLocation> clinicalTrialLocations) {
        List<ClinicalTrialLocation> clinicalTrialLocationsInterpretation = Lists.newArrayList();

        for (com.hartwig.hmftools.ckb.json.clinicaltrial.ClinicalTrialLocation location : clinicalTrialLocations) {
            clinicalTrialLocationsInterpretation.add(ImmutableClinicalTrialLocation.builder()
                    .nctId(location.nctId())
                    .facility(location.facility())
                    .city(location.city())
                    .country(location.country())
                    .status(location.status())
                    .state(location.state())
                    .zip(location.zip())
                    .contacts(extractClinicalTrialContacts(location.clinicalTrialContacts()))
                    .build());
        }
        return clinicalTrialLocationsInterpretation;
    }

    @NotNull
    private static List<com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ClinicalTrialContact> extractClinicalTrialContacts(
            @NotNull List<ClinicalTrialContact> contacts) {
        List<com.hartwig.hmftools.ckb.datamodelinterpretation.clinicaltrial.ClinicalTrialContact> clinicalTrialContacts =
                Lists.newArrayList();

        for (ClinicalTrialContact contact : contacts) {
            clinicalTrialContacts.add(ImmutableClinicalTrialContact.builder()
                    .name(contact.name())
                    .email(contact.email())
                    .phone(contact.phone())
                    .phoneExt(contact.phoneExt())
                    .role(contact.role())
                    .build());
        }
        return clinicalTrialContacts;
    }
}