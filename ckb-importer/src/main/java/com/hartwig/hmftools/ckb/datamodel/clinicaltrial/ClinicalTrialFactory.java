package com.hartwig.hmftools.ckb.datamodel.clinicaltrial;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.indication.IndicationFactory;
import com.hartwig.hmftools.ckb.datamodel.therapy.TherapyFactory;
import com.hartwig.hmftools.ckb.json.CkbJsonDatabase;
import com.hartwig.hmftools.ckb.json.clinicaltrial.JsonClinicalTrial;
import com.hartwig.hmftools.ckb.json.clinicaltrial.JsonContact;
import com.hartwig.hmftools.ckb.json.clinicaltrial.JsonLocation;
import com.hartwig.hmftools.ckb.json.clinicaltrial.JsonVariantRequirementDetail;
import com.hartwig.hmftools.ckb.json.common.ClinicalTrialInfo;

import org.jetbrains.annotations.NotNull;

public final class ClinicalTrialFactory {

    private ClinicalTrialFactory() {
    }

    @NotNull
    public static List<ClinicalTrial> extractClinicalTrials(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull List<ClinicalTrialInfo> clinicalTrialInfos) {
        List<ClinicalTrial> clinicalTrials = Lists.newArrayList();
        for (ClinicalTrialInfo clinicalTrialInfo : clinicalTrialInfos) {
            clinicalTrials.add(resolveClinicalTrial(ckbJsonDatabase, clinicalTrialInfo));
        }
        return clinicalTrials;
    }

    @NotNull
    private static ClinicalTrial resolveClinicalTrial(@NotNull CkbJsonDatabase ckbJsonDatabase,
            @NotNull ClinicalTrialInfo clinicalTrialInfo) {
        for (JsonClinicalTrial clinicalTrial : ckbJsonDatabase.clinicalTrials()) {
            if (clinicalTrialInfo.nctId().equals(clinicalTrial.nctId())) {
                return ImmutableClinicalTrial.builder()
                        .updateDate(clinicalTrial.updateDate())
                        .nctId(clinicalTrial.nctId())
                        .title(clinicalTrial.title())
                        .therapies(TherapyFactory.extractTherapies(ckbJsonDatabase, clinicalTrial.therapies()))
                        .indications(IndicationFactory.extractIndications(ckbJsonDatabase, clinicalTrial.indications()))
                        .phase(clinicalTrial.phase())
                        .recruitment(clinicalTrial.recruitment())
                        .ageGroups(clinicalTrial.ageGroups())
                        .gender(clinicalTrial.gender())
                        .variantRequirement(clinicalTrial.variantRequirements())
                        .sponsor(clinicalTrial.sponsors())
                        .variantRequirementDetails(convertRequirementDetails(clinicalTrial.variantRequirementDetails()))
                        .locations(convertLocations(clinicalTrial.locations()))
                        .build();
            }
        }

        throw new IllegalStateException("Could not resolve CKB clinical trial with nct '" + clinicalTrialInfo.nctId() + "'");
    }

    @NotNull
    private static List<VariantRequirementDetail> convertRequirementDetails(
            @NotNull List<JsonVariantRequirementDetail> jsonRequirementDetails) {
        List<VariantRequirementDetail> requirementDetails = Lists.newArrayList();
        for (JsonVariantRequirementDetail requirementDetail : jsonRequirementDetails) {
            requirementDetails.add(ImmutableVariantRequirementDetail.builder()
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

    @NotNull
    private static List<Location> convertLocations(@NotNull List<JsonLocation> jsonLocations) {
        List<Location> locations = Lists.newArrayList();
        for (JsonLocation location : jsonLocations) {
            locations.add(ImmutableLocation.builder()
                    .nctId(location.nctId())
                    .facility(location.facility())
                    .city(location.city())
                    .country(location.country())
                    .status(location.status())
                    .state(location.state())
                    .zip(location.zip())
                    .contacts(convertContacts(location.contacts()))
                    .build());
        }
        return locations;
    }

    @NotNull
    private static List<Contact> convertContacts(@NotNull List<JsonContact> jsonContacts) {
        List<Contact> contacts = Lists.newArrayList();
        for (JsonContact contact : jsonContacts) {
            contacts.add(ImmutableContact.builder()
                    .name(contact.name())
                    .email(contact.email())
                    .phone(contact.phone())
                    .phoneExt(contact.phoneExt())
                    .role(contact.role())
                    .build());
        }
        return contacts;
    }
}