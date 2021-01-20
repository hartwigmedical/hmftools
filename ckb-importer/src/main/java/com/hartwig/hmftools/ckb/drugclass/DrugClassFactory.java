package com.hartwig.hmftools.ckb.drugclass;

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
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class DrugClassFactory {

    private static final Logger LOGGER = LogManager.getLogger(DrugClassFactory.class);

    private DrugClassFactory() {
    }

    @NotNull
    public static List<DrugClass> readingDrugClasses(@NotNull String drugClassesDir) throws IOException {
        LOGGER.info("Start reading drug classes");

        List<DrugClass> drugClasses = Lists.newArrayList();
        File[] filesDrugClasses = new File(drugClassesDir).listFiles();

        if (filesDrugClasses != null) {
            LOGGER.info("The total files in the drug classes dir is {}", filesDrugClasses.length);

            for (File drugClass : filesDrugClasses) {
                JsonParser parser = new JsonParser();
                JsonReader reader = new JsonReader(new FileReader(drugClass));
                reader.setLenient(true);

                while (reader.peek() != JsonToken.END_DOCUMENT) {
                    JsonObject drugClassEntryObject = parser.parse(reader).getAsJsonObject();
                    JsonDatamodelChecker drugsClassChecker = DrugClassDataModelChecker.drugClassObjectChecker();
                    drugsClassChecker.check(drugClassEntryObject);

                    drugClasses.add(ImmutableDrugClass.builder()
                            .id(JsonFunctions.string(drugClassEntryObject, "id"))
                            .drugClass(JsonFunctions.string(drugClassEntryObject, "drugClass"))
                            .createDate(JsonFunctions.string(drugClassEntryObject, "createDate"))
                            .drugs(retrieveDrugs(drugClassEntryObject.getAsJsonArray("drugs")))
                            .treatmentApproaches(retrieveTreatmentApproaches(drugClassEntryObject.getAsJsonArray("treatmentApproaches")))
                            .build());
                }
            }
        }
        LOGGER.info("Finished reading drug classes");
        return drugClasses;
    }

    private static List<DrugClassDrug> retrieveDrugs(@NotNull JsonArray jsonArray) {
        List<DrugClassDrug> drugs = Lists.newArrayList();
        JsonDatamodelChecker drugsClassDrugChecker = DrugClassDataModelChecker.drugClassDrugsObjectChecker();
        for (JsonElement drug : jsonArray) {
            JsonObject drugsObject = drug.getAsJsonObject();
            drugsClassDrugChecker.check(drugsObject);

            drugs.add(ImmutableDrugClassDrug.builder()
                    .id(JsonFunctions.string(drugsObject, "id"))
                    .drugName(JsonFunctions.string(drugsObject, "drugName"))
                    .drugsTerms(JsonFunctions.optionalStringList(drugsObject, "terms"))
                    .build());
        }
        return drugs;
    }

    private static List<DrugClassTreatmentApproach> retrieveTreatmentApproaches(@NotNull JsonArray jsonArray) {
        List<DrugClassTreatmentApproach> treatmentApproaches = Lists.newArrayList();
        JsonDatamodelChecker drugsClassTreatmentApprochChecker = DrugClassDataModelChecker.drugClassTreatmentApproachesObjectChecker();

        for (JsonElement treatmentApproach : jsonArray) {
            JsonObject treatmentApproachObject = treatmentApproach.getAsJsonObject();
            drugsClassTreatmentApprochChecker.check(treatmentApproachObject);

            treatmentApproaches.add(ImmutableDrugClassTreatmentApproach.builder()
                    .id(JsonFunctions.string(treatmentApproachObject, "id"))
                    .name(JsonFunctions.string(treatmentApproachObject, "name"))
                    .profileName(JsonFunctions.string(treatmentApproachObject, "profileName"))
                    .build());
        }
        return treatmentApproaches;
    }
}
