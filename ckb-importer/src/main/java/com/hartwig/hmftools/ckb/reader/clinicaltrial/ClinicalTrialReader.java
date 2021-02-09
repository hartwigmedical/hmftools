package com.hartwig.hmftools.ckb.reader.clinicaltrial;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrial;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrialContact;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrialLocation;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ClinicalTrialVariantRequirementDetail;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ImmutableClinicalTrial;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ImmutableClinicalTrialContact;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ImmutableClinicalTrialLocation;
import com.hartwig.hmftools.ckb.datamodel.clinicaltrial.ImmutableClinicalTrialVariantRequirementDetail;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableIndicationInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableMolecularProfileInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableTherapyInfo;
import com.hartwig.hmftools.ckb.datamodel.common.IndicationInfo;
import com.hartwig.hmftools.ckb.datamodel.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.datamodel.common.TherapyInfo;
import com.hartwig.hmftools.ckb.reader.CkbJsonDirectoryReader;
import com.hartwig.hmftools.ckb.util.DateConverter;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ClinicalTrialReader extends CkbJsonDirectoryReader<ClinicalTrial> {

    public ClinicalTrialReader(@Nullable final Integer maxFilesToRead) {
        super(maxFilesToRead);
    }

    @NotNull
    @Override
    protected ClinicalTrial read(@NotNull final JsonObject object) {
        JsonDatamodelChecker clinicalTrialChecker = ClinicalTrialDataModelChecker.clinicalTrialObjectChecker();
        clinicalTrialChecker.check(object);

        return ImmutableClinicalTrial.builder()
                .nctId(JsonFunctions.string(object, "nctId"))
                .title(JsonFunctions.string(object, "title"))
                .phase(JsonFunctions.string(object, "phase"))
                .recruitment(JsonFunctions.string(object, "recruitment"))
                .therapy(extractClinicalTrialsTherapies(object.getAsJsonArray("therapies")))
                .ageGroup(JsonFunctions.stringList(object, "ageGroups"))
                .gender(JsonFunctions.optionalNullableString(object, "gender"))
                .variantRequirement(JsonFunctions.string(object, "variantRequirements"))
                .sponsors(JsonFunctions.optionalNullableString(object, "sponsors"))
                .updateDate(DateConverter.toDate(JsonFunctions.string(object, "updateDate")))
                .indication(extractClinicalTrialsIndications(object.getAsJsonArray("indications")))
                .variantRequirementDetail(extractClinicalTrialsVariantRequirementDetails(object.getAsJsonArray("variantRequirementDetails")))
                .clinicalTrialLocation(extractClinicalTrialsLocations(object.getAsJsonArray("clinicalTrialLocations")))
                .build();
    }

    @NotNull
    private static List<ClinicalTrialVariantRequirementDetail> extractClinicalTrialsVariantRequirementDetails(
            @NotNull JsonArray jsonArray) {
        List<ClinicalTrialVariantRequirementDetail> variantRequirementDetails = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrailVariantRequirementDetailChecker =
                ClinicalTrialDataModelChecker.clinicalTrialVariantRequirementDetailsObjectChecker();

        for (JsonElement variantRequirementDetail : jsonArray) {
            JsonObject variantRequirementDetailObject = variantRequirementDetail.getAsJsonObject();
            clinicalTrailVariantRequirementDetailChecker.check(variantRequirementDetailObject);

            variantRequirementDetails.add(ImmutableClinicalTrialVariantRequirementDetail.builder()
                    .molecularProfile(extractClinicalTrialsMolecularProfile(variantRequirementDetailObject.getAsJsonObject(
                            "molecularProfile")))
                    .requirementType(JsonFunctions.string(variantRequirementDetailObject, "requirementType"))
                    .build());
        }
        return variantRequirementDetails;
    }

    @NotNull
    private static MolecularProfileInfo extractClinicalTrialsMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker clinicalTrailMolecularProfileChecker =
                ClinicalTrialDataModelChecker.clinicalTrialMolecularProfileObjectChecker();
        clinicalTrailMolecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

    @NotNull
    private static List<TherapyInfo> extractClinicalTrialsTherapies(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> therapies = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrailTherapiesChecker = ClinicalTrialDataModelChecker.clinicalTrialTherapiesObjectChecker();

        for (JsonElement therapy : jsonArray) {
            JsonObject therapyObject = therapy.getAsJsonObject();
            clinicalTrailTherapiesChecker.check(therapyObject);

            therapies.add(ImmutableTherapyInfo.builder()
                    .id(JsonFunctions.integer(therapyObject, "id"))
                    .therapyName(JsonFunctions.string(therapyObject, "therapyName"))
                    .synonyms(JsonFunctions.optionalNullableString(therapyObject, "synonyms"))
                    .build());
        }
        return therapies;
    }

    @NotNull
    private static List<IndicationInfo> extractClinicalTrialsIndications(@NotNull JsonArray jsonArray) {
        List<IndicationInfo> clinicalTrialsIndications = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrailIndicationsChecker = ClinicalTrialDataModelChecker.clinicalTrialIndicationsObjectChecker();

        for (JsonElement clinicalTrialIndication : jsonArray) {
            JsonObject clinicalTrialIndicationObject = clinicalTrialIndication.getAsJsonObject();
            clinicalTrailIndicationsChecker.check(clinicalTrialIndicationObject);

            clinicalTrialsIndications.add(ImmutableIndicationInfo.builder()
                    .id(JsonFunctions.integer(clinicalTrialIndicationObject, "id"))
                    .name(JsonFunctions.string(clinicalTrialIndicationObject, "name"))
                    .source(JsonFunctions.string(clinicalTrialIndicationObject, "source"))
                    .build());
        }
        return clinicalTrialsIndications;
    }

    @NotNull
    private static List<ClinicalTrialLocation> extractClinicalTrialsLocations(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialLocation> clinicalTrialLocations = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrailLocationsChecker = ClinicalTrialDataModelChecker.clinicalTrialLocationsObjectChecker();

        for (JsonElement clinicalTrialLocation : jsonArray) {
            JsonObject clinicalTrailLocationObject = clinicalTrialLocation.getAsJsonObject();
            clinicalTrailLocationsChecker.check(clinicalTrailLocationObject);

            clinicalTrialLocations.add(ImmutableClinicalTrialLocation.builder()
                    .nctId(JsonFunctions.string(clinicalTrailLocationObject, "nctId"))
                    .facility(JsonFunctions.optionalNullableString(clinicalTrailLocationObject, "facility"))
                    .city(JsonFunctions.string(clinicalTrailLocationObject, "city"))
                    .country(JsonFunctions.string(clinicalTrailLocationObject, "country"))
                    .status(JsonFunctions.optionalNullableString(clinicalTrailLocationObject, "status"))
                    .state(JsonFunctions.optionalNullableString(clinicalTrailLocationObject, "state"))
                    .zip(JsonFunctions.optionalNullableString(clinicalTrailLocationObject, "zip"))
                    .clinicalTrialContact(extractClinicalTrialsContact(clinicalTrailLocationObject.getAsJsonArray("clinicalTrialContacts")))
                    .build());
        }
        return clinicalTrialLocations;
    }

    @NotNull
    private static List<ClinicalTrialContact> extractClinicalTrialsContact(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialContact> clinicalTrialContacts = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrailContactChecker = ClinicalTrialDataModelChecker.clinicalTrialContactObjectChecker();

        for (JsonElement clinicalTrialContact : jsonArray) {
            JsonObject clinicalTrialContactObject = clinicalTrialContact.getAsJsonObject();

            clinicalTrailContactChecker.check(clinicalTrialContactObject);
            clinicalTrialContacts.add(ImmutableClinicalTrialContact.builder()
                    .name(JsonFunctions.optionalNullableString(clinicalTrialContactObject, "name"))
                    .email(JsonFunctions.optionalNullableString(clinicalTrialContactObject, "email"))
                    .phone(JsonFunctions.optionalNullableString(clinicalTrialContactObject, "phone"))
                    .phoneExt(JsonFunctions.optionalNullableString(clinicalTrialContactObject, "phoneExt"))
                    .role(JsonFunctions.string(clinicalTrialContactObject, "role"))
                    .build());

        }
        return clinicalTrialContacts;
    }
}
