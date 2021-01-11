package com.hartwig.hmftools.ckb.clinicaltrial;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.hartwig.hmftools.ckb.CkbDataModelChecker;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ClinicalTrialFactory {

    private static final Logger LOGGER = LogManager.getLogger(ClinicalTrialFactory.class);

    private ClinicalTrialFactory() {

    }

    @NotNull
    public static List<ClinicalTrial> readingClinicalTrial(@NotNull String clinicalTrialDir) throws IOException {
        LOGGER.info("Start reading clinical trials");

        List<ClinicalTrial> clinicalTrials = Lists.newArrayList();
        File[] filesClinicalTrials = new File(clinicalTrialDir).listFiles();

        if (filesClinicalTrials != null) {
            LOGGER.info("The total files in the clinical trial dir is {}", filesClinicalTrials.length);

            for (File clinicalTrial : filesClinicalTrials) {
                JsonParser parser = new JsonParser();
                JsonReader reader = new JsonReader(new FileReader(clinicalTrial));
                reader.setLenient(true);

                while (reader.peek() != JsonToken.END_DOCUMENT) {
                    JsonObject clinicalTrialsEntryObject = parser.parse(reader).getAsJsonObject();
                    JsonDatamodelChecker clinicalTrailChecker = CkbDataModelChecker.clinicalTrialObjectChecker();
                    clinicalTrailChecker.check(clinicalTrialsEntryObject);

                    clinicalTrials.add(ImmutableClinicalTrial.builder()
                            .nctId(JsonFunctions.string(clinicalTrialsEntryObject, "nctId"))
                            .title(JsonFunctions.string(clinicalTrialsEntryObject, "title"))
                            .phase(JsonFunctions.string(clinicalTrialsEntryObject, "phase"))
                            .recruitment(JsonFunctions.string(clinicalTrialsEntryObject, "recruitment"))
                            .therapies(retrieveClinicalTrialsTherapies(clinicalTrialsEntryObject.getAsJsonArray("therapies")))
                            .ageGroups(JsonFunctions.stringList(clinicalTrialsEntryObject, "ageGroups"))
                            .gender(JsonFunctions.optionalNullableString(clinicalTrialsEntryObject, "gender"))
                            .variantRequirements(JsonFunctions.string(clinicalTrialsEntryObject, "variantRequirements"))
                            .sponsors(JsonFunctions.optionalNullableString(clinicalTrialsEntryObject, "sponsors"))
                            .updateDate(JsonFunctions.string(clinicalTrialsEntryObject, "updateDate"))
                            .indications(retrieveClinicalTrialsIndications(clinicalTrialsEntryObject.getAsJsonArray("indications")))
                            .variantRequirementDetails(retrieveClinicalTrialsVariantRequirementDetails(clinicalTrialsEntryObject.getAsJsonArray(
                                    "variantRequirementDetails")))
                            .clinicalTrialLocations(retrieveClinicalTrialsLocations(clinicalTrialsEntryObject.getAsJsonArray(
                                    "clinicalTrialLocations")))
                            .build());
                }
            }
        }
        LOGGER.info("Finished reading clinical trials");

        return clinicalTrials;
    }

    private static List<ClinicalTrialVariantRequirementDetails> retrieveClinicalTrialsVariantRequirementDetails(
            @NotNull JsonArray jsonArray) {
        List<ClinicalTrialVariantRequirementDetails> variantRequirementDetails = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrailVariantRequirementDetailChecker =
                CkbDataModelChecker.clinicalTrialVariantRequirementDetailsObjectChecker();
        for (JsonElement variantDetail : jsonArray) {
            clinicalTrailVariantRequirementDetailChecker.check(variantDetail.getAsJsonObject());
            variantRequirementDetails.add(ImmutableClinicalTrialVariantRequirementDetails.builder()
                    .molecularProfile(retrieveClinicalTrialsMolecularProfile(variantDetail.getAsJsonObject()
                            .getAsJsonObject("molecularProfile")))
                    .requirementType(JsonFunctions.string(variantDetail.getAsJsonObject(), "requirementType"))
                    .build());
        }
        return variantRequirementDetails;
    }

    private static ClinicalTrailMolecularProfile retrieveClinicalTrialsMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker clinicalTrailMolecularProfileChecker = CkbDataModelChecker.clinicalTrialMolecularProfileObjectChecker();
        clinicalTrailMolecularProfileChecker.check(jsonObject);
        return ImmutableClinicalTrailMolecularProfile.builder()
                .id(JsonFunctions.nullableString(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();

    }

    private static List<Therapies> retrieveClinicalTrialsTherapies(@NotNull JsonArray jsonArray) {
        List<Therapies> therapies = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrailTherapiesChecker = CkbDataModelChecker.clinicalTrialTherapiesObjectChecker();
        for (JsonElement therapy : jsonArray) {
            clinicalTrailTherapiesChecker.check(therapy.getAsJsonObject());

            therapies.add(ImmutableTherapies.builder()
                    .id(JsonFunctions.string(therapy.getAsJsonObject(), "id"))
                    .therapyName(JsonFunctions.string(therapy.getAsJsonObject(), "therapyName"))
                    .synonyms(JsonFunctions.optionalNullableString(therapy.getAsJsonObject(), "synonyms"))
                    .build());
        }
        return therapies;
    }

    private static List<Indications> retrieveClinicalTrialsIndications(@NotNull JsonArray jsonArray) {
        List<Indications> indications = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrailIndicationsChecker = CkbDataModelChecker.clinicalTrialIndicationsObjectChecker();

        for (JsonElement indication : jsonArray) {
            clinicalTrailIndicationsChecker.check(indication.getAsJsonObject());
            indications.add(ImmutableIndications.builder()
                    .id(JsonFunctions.string(indication.getAsJsonObject(), "id"))
                    .name(JsonFunctions.string(indication.getAsJsonObject(), "name"))
                    .source(JsonFunctions.string(indication.getAsJsonObject(), "source"))
                    .build());
        }
        return indications;
    }

    private static List<ClinicalTrialLocations> retrieveClinicalTrialsLocations(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialLocations> locations = Lists.newArrayList();
        JsonDatamodelChecker clinicalTrailLocationsChecker = CkbDataModelChecker.clinicalTrialLocationsObjectChecker();

        for (JsonElement location : jsonArray) {
            clinicalTrailLocationsChecker.check(location.getAsJsonObject());
            locations.add(ImmutableClinicalTrialLocations.builder()
                    .nctId(JsonFunctions.string(location.getAsJsonObject(), "nctId"))
                    .facility(JsonFunctions.optionalNullableString(location.getAsJsonObject(), "facility"))
                    .city(JsonFunctions.string(location.getAsJsonObject(), "city"))
                    .country(JsonFunctions.string(location.getAsJsonObject(), "country"))
                    .status(JsonFunctions.optionalNullableString(location.getAsJsonObject(), "status"))
                    .state(JsonFunctions.optionalNullableString(location.getAsJsonObject(), "state"))
                    .zip(JsonFunctions.optionalNullableString(location.getAsJsonObject(), "zip"))
                    .clinicalTrialContacts(retrieveClinicalTrialsContact(location.getAsJsonObject()))
                    .build());
        }
        return locations;
    }

    private static List<ClinicalTrialContacts> retrieveClinicalTrialsContact(@NotNull JsonObject jsonObject) {
        List<ClinicalTrialContacts> contacts = Lists.newArrayList();

        JsonDatamodelChecker clinicalTrailContactChecker = CkbDataModelChecker.clinicalTrialContactObjectChecker();
        JsonArray arrayContact = jsonObject.getAsJsonArray("clinicalTrialContacts");

        for (JsonElement contact : arrayContact) {
            clinicalTrailContactChecker.check(contact.getAsJsonObject());
            contacts.add(ImmutableClinicalTrialContacts.builder()
                    .name(JsonFunctions.optionalNullableString(contact.getAsJsonObject(), "name"))
                    .email(JsonFunctions.optionalNullableString(contact.getAsJsonObject(), "email"))
                    .phone(JsonFunctions.optionalNullableString(contact.getAsJsonObject(), "phone"))
                    .phoneExt(JsonFunctions.optionalNullableString(contact.getAsJsonObject(), "phoneExt"))
                    .role(JsonFunctions.string(contact.getAsJsonObject(), "role"))
                    .build());

        }
        return contacts;

    }
}
