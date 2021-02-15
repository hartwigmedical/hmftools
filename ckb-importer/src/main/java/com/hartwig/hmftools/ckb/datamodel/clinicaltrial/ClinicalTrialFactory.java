package com.hartwig.hmftools.ckb.datamodel.clinicaltrial;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.ImmutableCkbEntry;
import com.hartwig.hmftools.ckb.datamodel.common.CommonInterpretationFactory;
import com.hartwig.hmftools.ckb.datamodel.common.molecularprofile.MolecularProfileInterpretationFactory;
import com.hartwig.hmftools.ckb.datamodel.common.therapy.TherapyInterpretationFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.clinicaltrial.JsonClinicalTrial;
import com.hartwig.hmftools.ckb.json.clinicaltrial.JsonClinicalTrialContact;
import com.hartwig.hmftools.ckb.json.clinicaltrial.JsonClinicalTrialLocation;
import com.hartwig.hmftools.ckb.json.common.ClinicalTrialInfo;
import com.hartwig.hmftools.ckb.json.common.IndicationInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
import com.hartwig.hmftools.ckb.json.molecularprofile.JsonMolecularProfile;
import com.hartwig.hmftools.ckb.json.therapy.JsonTherapy;

import org.jetbrains.annotations.NotNull;

public final class ClinicalTrialFactory {

    private ClinicalTrialFactory() {
    }

    public static void interpretClinicalTrials(@NotNull JsonMolecularProfile molecularProfile, @NotNull CkbJsonDatabase ckbEntry,
            @NotNull ImmutableCkbEntry.Builder outputBuilder) {
        for (ClinicalTrialInfo clinicalTrialInfo : molecularProfile.variantAssociatedClinicalTrials()) {
            ImmutableClinicalTrial.Builder outputBuilderClinical = ImmutableClinicalTrial.builder();

            for (JsonClinicalTrial clinicalTrial : ckbEntry.clinicalTrials()) {
                if (clinicalTrialInfo.nctId().equals(clinicalTrial.nctId())) {
                    outputBuilderClinical.nctId(clinicalTrial.nctId())
                            .title(clinicalTrial.title())
                            .phase(clinicalTrial.phase())
                            .recruitment(clinicalTrial.recruitment())
                            .ageGroups(clinicalTrial.ageGroups())
                            .gender(clinicalTrial.gender())
                            .variantRequirement(clinicalTrial.variantRequirements())
                            .sponsor(clinicalTrial.sponsors())
                            .updateDate(clinicalTrial.updateDate())
                            .clinicalTrialVariantRequirementDetails(MolecularProfileInterpretationFactory.extractProfileNameClinicalTrial(
                                    clinicalTrial.variantRequirementDetails(),
                                    molecularProfile,
                                    ckbEntry))
                            .locations(extractClinicalTrialLocations(clinicalTrial.clinicalTrialLocations()));

                    for (IndicationInfo indicationInfo : clinicalTrial.indications()) {
                        outputBuilderClinical.addIndications(CommonInterpretationFactory.extractIndication(ckbEntry, indicationInfo));
                    }

                    for (TherapyInfo therapyInfo : clinicalTrial.therapies()) {
                        for (JsonTherapy therapy : ckbEntry.therapies()) {
                            if (therapyInfo.id() == therapy.id()) {
                                outputBuilderClinical.addTherapies(TherapyInterpretationFactory.extractTherapy(therapy,
                                        ckbEntry,
                                        molecularProfile));
                            }
                        }
                    }
                }
            }
            outputBuilder.addClinicalTrials(outputBuilderClinical.build());
        }
    }

    @NotNull
    private static List<ClinicalTrialLocation> extractClinicalTrialLocations(
            @NotNull List<JsonClinicalTrialLocation> clinicalTrialLocations) {
        List<ClinicalTrialLocation> clinicalTrialLocationsInterpretation = Lists.newArrayList();

        for (JsonClinicalTrialLocation location : clinicalTrialLocations) {
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
    private static List<com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrialContact> extractClinicalTrialContacts(
            @NotNull List<JsonClinicalTrialContact> contacts) {
        List<com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrialContact> clinicalTrialContacts = Lists.newArrayList();

        for (JsonClinicalTrialContact contact : contacts) {
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