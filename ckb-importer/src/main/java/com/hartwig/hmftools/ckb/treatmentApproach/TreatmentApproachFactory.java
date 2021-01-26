package com.hartwig.hmftools.ckb.treatmentApproach;

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
import com.hartwig.hmftools.ckb.common.ImmutableReferenceInfo;
import com.hartwig.hmftools.ckb.common.ImmutableTherapyInfo;
import com.hartwig.hmftools.ckb.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.common.TherapyInfo;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class TreatmentApproachFactory {

    private static final Logger LOGGER = LogManager.getLogger(TreatmentApproachFactory.class);

    private TreatmentApproachFactory() {

    }

    @NotNull
    public static List<TreatmentApproach> readingTreatmentApproch(@NotNull String treatmentApprochDir) throws IOException {
        LOGGER.info("Start reading treatment approach");

        List<TreatmentApproach> treatmentAppraoch = Lists.newArrayList();
        File[] filesTreatmentApproch = new File(treatmentApprochDir).listFiles();

        if (filesTreatmentApproch != null) {
            LOGGER.info("The total files in the treatment approch dir is {}", filesTreatmentApproch.length);

            for (File treatmentApproch : filesTreatmentApproch) {
                JsonParser parser = new JsonParser();
                JsonReader reader = new JsonReader(new FileReader(treatmentApproch));
                reader.setLenient(true);

                while (reader.peek() != JsonToken.END_DOCUMENT) {
                    JsonObject treatmentApprochEntryObject = parser.parse(reader).getAsJsonObject();
                    JsonDatamodelChecker treatmentApprochObjectChecker = TreatmentApprochDataModelChecker.treatmentApprochObjectChecker();
                    treatmentApprochObjectChecker.check(treatmentApprochEntryObject);

                    treatmentAppraoch.add(ImmutableTreatmentApproach.builder()
                            .id(JsonFunctions.string(treatmentApprochEntryObject, "id"))
                            .name(JsonFunctions.string(treatmentApprochEntryObject, "name"))
                            .profileName(JsonFunctions.string(treatmentApprochEntryObject, "profileName"))
                            .drugClass(treatmentApprochEntryObject.has("drugClass") && !treatmentApprochEntryObject.get("drugClass")
                                    .isJsonNull() ? createDrugClass(treatmentApprochEntryObject.getAsJsonObject("drugClass")) : null)
                            .therapy(treatmentApprochEntryObject.has("therapy") && !treatmentApprochEntryObject.get("therapy").isJsonNull()
                                    ? extractTherapy(treatmentApprochEntryObject.getAsJsonObject("therapy"))
                                    : null)
                            .reference(extractReference(treatmentApprochEntryObject.getAsJsonArray("references")))
                            .createDate(JsonFunctions.string(treatmentApprochEntryObject, "createDate"))
                            .updateDate(JsonFunctions.string(treatmentApprochEntryObject, "updateDate"))
                            .build());
                }
                reader.close();
            }
        }
        LOGGER.info("Finished reading treatment approach");

        return treatmentAppraoch;
    }

    @NotNull
    public static TreatmentApprochDrugClass createDrugClass(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker drugClassObjectChecker = TreatmentApprochDataModelChecker.drugClassObjectChecker();
        drugClassObjectChecker.check(jsonObject);

        return ImmutableTreatmentApprochDrugClass.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .drugClass(JsonFunctions.string(jsonObject, "drugClass"))
                .build();
    }

    @NotNull
    public static TherapyInfo extractTherapy(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker therapyObjectChecker = TreatmentApprochDataModelChecker.therapyObjectChecker();
        therapyObjectChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.nullableString(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    public static List<ReferenceInfo> extractReference(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> references = Lists.newArrayList();
        JsonDatamodelChecker referenceObjectChecker = TreatmentApprochDataModelChecker.referenceObjectChecker();

        for (JsonElement reference : jsonArray) {
            JsonObject referenceJsonObject = reference.getAsJsonObject();
            referenceObjectChecker.check(referenceJsonObject);

            references.add(ImmutableReferenceInfo.builder()
                    .id(JsonFunctions.string(referenceJsonObject, "id"))
                    .pubMedId(JsonFunctions.nullableString(referenceJsonObject, "pubMedId"))
                    .title(JsonFunctions.string(referenceJsonObject, "title"))
                    .url(JsonFunctions.nullableString(referenceJsonObject, "url"))
                    .build());
        }

        return references;
    }
}
