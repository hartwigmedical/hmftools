package com.hartwig.hmftools.ckb.datamodel.clinicaltrial;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.common.CommonInterpretationFactory;
import com.hartwig.hmftools.ckb.datamodel.common.molecularprofile.MolecularProfileInterpretationFactory;
import com.hartwig.hmftools.ckb.datamodel.common.therapy.TherapyFactory;
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

    @NotNull
    public static List<ClinicalTrial> interpretClinicalTrials(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull JsonMolecularProfile molecularProfile) {
        List<ClinicalTrial> clinicalTrials = Lists.newArrayList();
        for (ClinicalTrialInfo clinicalTrialInfo : molecularProfile.variantAssociatedClinicalTrials()) {
            ImmutableClinicalTrial.Builder outputBuilderClinical = ImmutableClinicalTrial.builder();

            for (JsonClinicalTrial clinicalTrial : ckbJsonDatabase.clinicalTrials()) {
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
                                    ckbJsonDatabase,
                                    clinicalTrial.variantRequirementDetails(),
                                    molecularProfile))
                            .locations(extractClinicalTrialLocations(clinicalTrial.clinicalTrialLocations()));

                    for (IndicationInfo indicationInfo : clinicalTrial.indications()) {
                        outputBuilderClinical.addIndications(CommonInterpretationFactory.extractIndication(ckbJsonDatabase,
                                indicationInfo));
                    }

                    for (TherapyInfo therapyInfo : clinicalTrial.therapies()) {
                        for (JsonTherapy therapy : ckbJsonDatabase.therapies()) {
                            if (therapyInfo.id() == therapy.id()) {
                                outputBuilderClinical.addTherapies(TherapyFactory.extractTherapy(ckbJsonDatabase,
                                        therapy,
                                        molecularProfile));
                            }
                        }
                    }
                }
            }
            clinicalTrials.add(outputBuilderClinical.build());
        }

        return clinicalTrials;
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
    private static List<ClinicalTrialContact> extractClinicalTrialContacts(@NotNull List<JsonClinicalTrialContact> contacts) {
        List<ClinicalTrialContact> clinicalTrialContacts = Lists.newArrayList();

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