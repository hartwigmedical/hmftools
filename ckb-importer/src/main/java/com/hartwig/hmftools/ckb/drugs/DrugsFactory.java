package com.hartwig.hmftools.ckb.drugs;

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
import org.jetbrains.annotations.NotNull;

public class DrugsFactory {

    private static final Logger LOGGER = LogManager.getLogger(DrugsFactory.class);

    private DrugsFactory() {
    }

    @NotNull
    public static List<Drugs> readingDrugs(@NotNull String drugsDir) throws IOException {
        LOGGER.info("Start reading drugs");

        List<Drugs> drugs = Lists.newArrayList();
        File[] filesDrugs = new File(drugsDir).listFiles();

        if (filesDrugs != null) {
            LOGGER.info("The total files in the drugs dir is {}", filesDrugs.length);

            for (File drug : filesDrugs) {
                JsonParser parser = new JsonParser();
                JsonReader reader = new JsonReader(new FileReader(drug));
                reader.setLenient(true);

                while (reader.peek() != JsonToken.END_DOCUMENT) {
                    JsonObject drugsEntryObject = parser.parse(reader).getAsJsonObject();
                    drugs.add(ImmutableDrugs.builder()
                            .id(JsonFunctions.string(drugsEntryObject, "id"))
                            .drugName(JsonFunctions.string(drugsEntryObject, "drugName"))
                            .terms(JsonFunctions.stringList(drugsEntryObject, "terms"))
                            .synonyms(JsonFunctions.stringList(drugsEntryObject, "synonyms"))
                            .tradeName(JsonFunctions.nullableString(drugsEntryObject, "tradeName"))
                            .drugDescriptions(extractDrugDescriptions(drugsEntryObject.getAsJsonArray("drugDescriptions")))
                            .drugClasses(extractDrugsClasses(drugsEntryObject.getAsJsonArray("drugClasses")))
                            .casRegistryNum(JsonFunctions.nullableString(drugsEntryObject, "casRegistryNum"))
                            .nctId(JsonFunctions.nullableString(drugsEntryObject, "ncitId"))
                            .createDate(JsonFunctions.string(drugsEntryObject, "createDate"))
                            .clinicalTrials(Lists.newArrayList())
                            .evidence(Lists.newArrayList())
                            .therapies(Lists.newArrayList())
                            .globalApprovaStatus(Lists.newArrayList())
                            .build());
                }
            }
        }
        return drugs;
    }

    @NotNull
    public static List<DrugDescription> extractDrugDescriptions(@NotNull JsonArray jsonArray) {
        List<DrugDescription> drugDescriptions = Lists.newArrayList();

        for (JsonElement drugDescription : jsonArray) {
            drugDescriptions.add(ImmutableDrugDescription.builder()
                    .description(JsonFunctions.string(drugDescription.getAsJsonObject(), "description"))
                    .references(Lists.newArrayList())
                    .build());
        }
        return drugDescriptions;
    }

    @NotNull
    public static List<DrugClasses> extractDrugsClasses(@NotNull JsonArray jsonArray) {
        List<DrugClasses> drugClasses = Lists.newArrayList();
        for (JsonElement drugClass : jsonArray) {
            drugClasses.add(ImmutableDrugClasses.builder()
                    .id(JsonFunctions.string(drugClass.getAsJsonObject(), "id"))
                    .drugClass(JsonFunctions.string(drugClass.getAsJsonObject(), "drugClass"))
                    .build());
        }
        return drugClasses;
    }
}
