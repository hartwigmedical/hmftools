package com.hartwig.hmftools.ckb.datamodel.clinicaltrial;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.indication.IndicationFactory;
import com.hartwig.hmftools.ckb.datamodel.therapy.TherapyFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.clinicaltrial.JsonClinicalTrial;
import com.hartwig.hmftools.ckb.json.clinicaltrial.JsonClinicalTrialContact;
import com.hartwig.hmftools.ckb.json.clinicaltrial.JsonClinicalTrialLocation;
import com.hartwig.hmftools.ckb.json.clinicaltrial.JsonClinicalTrialVariantRequirementDetail;
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
                            .updateDate(clinicalTrial.updateDate())
                            .title(clinicalTrial.title())
                            .phase(clinicalTrial.phase())
                            .recruitment(clinicalTrial.recruitment())
                            .ageGroups(clinicalTrial.ageGroups())
                            .gender(clinicalTrial.gender())
                            .variantRequirement(clinicalTrial.variantRequirements())
                            .sponsor(clinicalTrial.sponsors())
                            .variantRequirementDetails(convertRequirementDetails(clinicalTrial.variantRequirementDetails()))
                            .locations(convertLocations(clinicalTrial.clinicalTrialLocations()));

                    for (IndicationInfo indicationInfo : clinicalTrial.indications()) {
                        outputBuilderClinical.addIndications(IndicationFactory.extractIndication(ckbJsonDatabase, indicationInfo));
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
    private static List<ClinicalTrialLocation> convertLocations(@NotNull List<JsonClinicalTrialLocation> jsonLocations) {
        List<ClinicalTrialLocation> locations = Lists.newArrayList();

        for (JsonClinicalTrialLocation location : jsonLocations) {
            locations.add(ImmutableClinicalTrialLocation.builder()
                    .nctId(location.nctId())
                    .facility(location.facility())
                    .city(location.city())
                    .country(location.country())
                    .status(location.status())
                    .state(location.state())
                    .zip(location.zip())
                    .contacts(convertContacts(location.clinicalTrialContacts()))
                    .build());
        }
        return locations;
    }

    @NotNull
    private static List<ClinicalTrialContact> convertContacts(@NotNull List<JsonClinicalTrialContact> jsonContacts) {
        List<ClinicalTrialContact> contacts = Lists.newArrayList();

        for (JsonClinicalTrialContact contact : jsonContacts) {
            contacts.add(ImmutableClinicalTrialContact.builder()
                    .name(contact.name())
                    .email(contact.email())
                    .phone(contact.phone())
                    .phoneExt(contact.phoneExt())
                    .role(contact.role())
                    .build());
        }
        return contacts;
    }

    @NotNull
    public static List<ClinicalTrialVariantRequirementDetail> convertRequirementDetails(
            @NotNull List<JsonClinicalTrialVariantRequirementDetail> jsonRequirementDetails) {
        List<ClinicalTrialVariantRequirementDetail> requirementDetails = Lists.newArrayList();
        for (JsonClinicalTrialVariantRequirementDetail requirementDetail : jsonRequirementDetails) {
            requirementDetails.add(ImmutableClinicalTrialVariantRequirementDetail.builder()
                    .profileId(requirementDetail.molecularProfile().id())
                    .requirementType(requirementDetail.requirementType())
                    .build());
            //            if (requirementDetail.requirementType().equals("excluded")) { // variant is excluded from enrollment
            //                requirementDetails.add(extractClinicalTrialVariantRequirementDetails(requirementDetail).build());
            //            }
            //
            //            if (requirementDetail.requirementType().equals("required")) { // variant is requirement for enrollment
            //                requirementDetails.add(extractClinicalTrialVariantRequirementDetails(requirementDetail).build());
            //            }
            //
            //            // variant is required or excluded for a subset of the enrollment population
            //            if (requirementDetail.requirementType().equals("partial - required")) {
            //                if (requirementDetail.molecularProfile().id() == molecularProfileDir.id()) {
            //                    requirementDetails.add(extractClinicalTrialVariantRequirementDetails(requirementDetail).build());
            //                }
            //            }
        }
        return requirementDetails;
    }
}