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
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public class ClinicalTrialFactory {

    private static final Logger LOGGER = LogManager.getLogger(ClinicalTrialFactory.class);

    private ClinicalTrialFactory(){

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
                    clinicalTrials.add(ImmutableClinicalTrial.builder()
                            .nctId(JsonFunctions.string(clinicalTrialsEntryObject, "nctId"))
                            .title(JsonFunctions.string(clinicalTrialsEntryObject, "title"))
                            .phase(JsonFunctions.string(clinicalTrialsEntryObject, "phase"))
                            .recruitment(JsonFunctions.string(clinicalTrialsEntryObject, "recruitment"))
                            .therapies(retrieveClinicalTrialsTherapies(
                                    clinicalTrialsEntryObject.getAsJsonArray("therapies")))
                            .ageGroups(Lists.newArrayList())
                            .gender(JsonFunctions.optionalNullableString(clinicalTrialsEntryObject, "gender"))
                            .variantRequirements(JsonFunctions.string(clinicalTrialsEntryObject, "variantRequirements"))
                            .sponsors(JsonFunctions.optionalNullableString(clinicalTrialsEntryObject, "sponsors"))
                            .updateDate(JsonFunctions.string(clinicalTrialsEntryObject, "updateDate"))
                            .indications(retrieveClinicalTrialsIndications(clinicalTrialsEntryObject.getAsJsonArray("indications")))
                            .variantRequirementDetails(Lists.newArrayList())
                            .clinicalTrialLocations(retrieveClinicalTrialsLocations(clinicalTrialsEntryObject.getAsJsonArray(
                                    "clinicalTrialLocations")))
                            .build());
                }
            }
        } LOGGER.info("Finished reading clinical trials");

        return clinicalTrials;
    }

    private static List<Therapies> retrieveClinicalTrialsTherapies(@NotNull JsonArray jsonArray) {
        List<Therapies> therapies = Lists.newArrayList();
        for (JsonElement therapy : jsonArray) {
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
        for (JsonElement indication : jsonArray) {
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

        for (JsonElement location : jsonArray) {
            locations.add(ImmutableClinicalTrialLocations.builder()
                    .nctId(JsonFunctions.string(location.getAsJsonObject(), "nctId"))
                    .facility(JsonFunctions.optionalNullableString(location.getAsJsonObject(), "facility"))
                    .city(JsonFunctions.string(location.getAsJsonObject(), "city"))
                    .country(JsonFunctions.string(location.getAsJsonObject(), "country"))
                    .status(JsonFunctions.optionalNullableString(location.getAsJsonObject(), "status"))
                    .state(JsonFunctions.optionalNullableString(location.getAsJsonObject(), "state"))
                    .zip(JsonFunctions.optionalNullableString(location.getAsJsonObject(), "zip"))
                    .clinicalTrialContacts(Lists.newArrayList())
                    .build());
        }
        return locations;
    }
}
