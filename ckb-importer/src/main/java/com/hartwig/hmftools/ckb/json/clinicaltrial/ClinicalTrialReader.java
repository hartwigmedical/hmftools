package com.hartwig.hmftools.ckb.json.clinicaltrial;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.ckb.json.CkbJsonDirectoryReader;
import com.hartwig.hmftools.ckb.json.common.ImmutableIndicationInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableMolecularProfileInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableTherapyInfo;
import com.hartwig.hmftools.ckb.json.common.IndicationInfo;
import com.hartwig.hmftools.ckb.json.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
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
                .therapies(extractTherapies(object.getAsJsonArray("therapies")))
                .ageGroups(JsonFunctions.stringList(object, "ageGroups"))
                .gender(JsonFunctions.optionalNullableString(object, "gender"))
                .variantRequirements(JsonFunctions.string(object, "variantRequirements"))
                .sponsors(JsonFunctions.optionalNullableString(object, "sponsors"))
                .updateDate(DateConverter.toDate(JsonFunctions.string(object, "updateDate")))
                .indications(extractIndications(object.getAsJsonArray("indications")))
                .variantRequirementDetails(extractVariantRequirementDetails(object.getAsJsonArray("variantRequirementDetails")))
                .clinicalTrialLocations(extractClinicalTrialsLocations(object.getAsJsonArray("clinicalTrialLocations")))
                .build();
    }

    @NotNull
    private static List<TherapyInfo> extractTherapies(@NotNull JsonArray jsonArray) {
        List<TherapyInfo> therapies = Lists.newArrayList();
        JsonDatamodelChecker therapyChecker = ClinicalTrialDataModelChecker.therapyObjectChecker();

        for (JsonElement therapy : jsonArray) {
            JsonObject therapyObject = therapy.getAsJsonObject();
            therapyChecker.check(therapyObject);

            therapies.add(ImmutableTherapyInfo.builder()
                    .id(JsonFunctions.integer(therapyObject, "id"))
                    .therapyName(JsonFunctions.string(therapyObject, "therapyName"))
                    .synonyms(JsonFunctions.optionalNullableString(therapyObject, "synonyms"))
                    .build());
        }
        return therapies;
    }

    @NotNull
    private static List<IndicationInfo> extractIndications(@NotNull JsonArray jsonArray) {
        List<IndicationInfo> indications = Lists.newArrayList();
        JsonDatamodelChecker indicationChecker = ClinicalTrialDataModelChecker.indicationObjectChecker();

        for (JsonElement indication : jsonArray) {
            JsonObject indicationObject = indication.getAsJsonObject();
            indicationChecker.check(indicationObject);

            indications.add(ImmutableIndicationInfo.builder()
                    .id(JsonFunctions.string(indicationObject, "id"))
                    .name(JsonFunctions.string(indicationObject, "name"))
                    .source(JsonFunctions.string(indicationObject, "source"))
                    .build());
        }
        return indications;
    }

    @NotNull
    private static MolecularProfileInfo extractMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker molecularProfileChecker = ClinicalTrialDataModelChecker.molecularProfileObjectChecker();
        molecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

    @NotNull
    private static List<ClinicalTrialVariantRequirementDetail> extractVariantRequirementDetails(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialVariantRequirementDetail> variantRequirementDetails = Lists.newArrayList();
        JsonDatamodelChecker variantRequirementDetailChecker =
                ClinicalTrialDataModelChecker.variantRequirementDetailObjectChecker();

        for (JsonElement variantRequirementDetail : jsonArray) {
            JsonObject variantRequirementDetailObject = variantRequirementDetail.getAsJsonObject();
            variantRequirementDetailChecker.check(variantRequirementDetailObject);

            variantRequirementDetails.add(ImmutableClinicalTrialVariantRequirementDetail.builder()
                    .molecularProfile(extractMolecularProfile(variantRequirementDetailObject.getAsJsonObject("molecularProfile")))
                    .requirementType(JsonFunctions.string(variantRequirementDetailObject, "requirementType"))
                    .build());
        }
        return variantRequirementDetails;
    }

    @NotNull
    private static List<ClinicalTrialLocation> extractClinicalTrialsLocations(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialLocation> locations = Lists.newArrayList();
        JsonDatamodelChecker locationChecker = ClinicalTrialDataModelChecker.locationObjectChecker();

        for (JsonElement location : jsonArray) {
            JsonObject locationObject = location.getAsJsonObject();
            locationChecker.check(locationObject);

            locations.add(ImmutableClinicalTrialLocation.builder()
                    .nctId(JsonFunctions.string(locationObject, "nctId"))
                    .facility(JsonFunctions.optionalNullableString(locationObject, "facility"))
                    .city(JsonFunctions.string(locationObject, "city"))
                    .country(JsonFunctions.string(locationObject, "country"))
                    .status(JsonFunctions.optionalNullableString(locationObject, "status"))
                    .state(JsonFunctions.optionalNullableString(locationObject, "state"))
                    .zip(JsonFunctions.optionalNullableString(locationObject, "zip"))
                    .clinicalTrialContacts(extractContacts(locationObject.getAsJsonArray("clinicalTrialContacts")))
                    .build());
        }
        return locations;
    }

    @NotNull
    private static List<ClinicalTrialContact> extractContacts(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialContact> contacts = Lists.newArrayList();
        JsonDatamodelChecker contactChecker = ClinicalTrialDataModelChecker.contactObjectChecker();

        for (JsonElement contact : jsonArray) {
            JsonObject contactObject = contact.getAsJsonObject();

            contactChecker.check(contactObject);
            contacts.add(ImmutableClinicalTrialContact.builder()
                    .name(JsonFunctions.optionalNullableString(contactObject, "name"))
                    .email(JsonFunctions.optionalNullableString(contactObject, "email"))
                    .phone(JsonFunctions.optionalNullableString(contactObject, "phone"))
                    .phoneExt(JsonFunctions.optionalNullableString(contactObject, "phoneExt"))
                    .role(JsonFunctions.string(contactObject, "role"))
                    .build());
        }
        return contacts;
    }
}
