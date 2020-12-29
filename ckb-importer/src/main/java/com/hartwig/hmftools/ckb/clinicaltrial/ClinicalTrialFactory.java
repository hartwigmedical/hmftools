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
                            .nctId(clinicalTrialsEntryObject.getAsJsonPrimitive("nctId").getAsString())
                            .title(clinicalTrialsEntryObject.getAsJsonPrimitive("title").getAsString())
                            .phase(clinicalTrialsEntryObject.getAsJsonPrimitive("phase").getAsString())
                            .recruitment(clinicalTrialsEntryObject.getAsJsonPrimitive("recruitment").getAsString())
                            .therapies(clinicalTrialsEntryObject.has("therapies") ? retrieveClinicalTrialsTherpaies(
                                    clinicalTrialsEntryObject.getAsJsonArray("therapies")) : null)
                            .ageGroups(Lists.newArrayList())
                            .gender(Strings.EMPTY)
                            .variantRequirements(clinicalTrialsEntryObject.getAsJsonPrimitive("variantRequirements").getAsString())
                            .sponsors(Strings.EMPTY)
                            .updateDate(clinicalTrialsEntryObject.getAsJsonPrimitive("updateDate").getAsString())
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

    private static List<Therapies> retrieveClinicalTrialsTherpaies(@NotNull JsonArray jsonArray) {
        List<Therapies> therapies = Lists.newArrayList();
        for (JsonElement therapy : jsonArray) {
            therapies.add(ImmutableTherapies.builder()
                    .id(therapy.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .therapyName(therapy.getAsJsonObject().getAsJsonPrimitive("therapyName").getAsString())
                    .synonyms(Strings.EMPTY)
                    .build());
        }
        return therapies;
    }

    private static List<Indications> retrieveClinicalTrialsIndications(@NotNull JsonArray jsonArray) {
        List<Indications> indications = Lists.newArrayList();
        for (JsonElement indication : jsonArray) {
            indications.add(ImmutableIndications.builder()
                    .id(indication.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .name(indication.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .source(indication.getAsJsonObject().getAsJsonPrimitive("source").getAsString())
                    .build());
        }
        return indications;
    }

    private static List<ClinicalTrialLocations> retrieveClinicalTrialsLocations(@NotNull JsonArray jsonArray) {
        List<ClinicalTrialLocations> locations = Lists.newArrayList();

        for (JsonElement location : jsonArray) {
            locations.add(ImmutableClinicalTrialLocations.builder()
                    .nctId(location.getAsJsonObject().getAsJsonPrimitive("nctId").getAsString())
                    .facility(Strings.EMPTY)
                    .city(location.getAsJsonObject().getAsJsonPrimitive("city").getAsString())
                    .country(location.getAsJsonObject().getAsJsonPrimitive("country").getAsString())
                    .status(Strings.EMPTY)
                    .state(Strings.EMPTY)
                    .zip(Strings.EMPTY)
                    .clinicalTrialContacts(Lists.newArrayList())

                    .build());
        }
        return locations;
    }

}
