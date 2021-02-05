package com.hartwig.hmftools.ckb.reader.treatmentapproch;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.hartwig.hmftools.ckb.datamodel.common.DrugClassInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableDrugClassInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableReferenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableTherapyInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.TherapyInfo;
import com.hartwig.hmftools.ckb.datamodel.treatmentapproach.ImmutableTreatmentApproach;
import com.hartwig.hmftools.ckb.datamodel.treatmentapproach.TreatmentApproach;
import com.hartwig.hmftools.ckb.util.DateConverter;
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
    public static List<TreatmentApproach> readingTreatmentApproch(@NotNull String treatmentApprochDir) throws IOException, ParseException {
        LOGGER.info("Start reading treatment approach dir");

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
                            .id(JsonFunctions.integer(treatmentApprochEntryObject, "id"))
                            .name(JsonFunctions.string(treatmentApprochEntryObject, "name"))
                            .profileName(JsonFunctions.string(treatmentApprochEntryObject, "profileName"))
                            .drugClass(treatmentApprochEntryObject.has("drugClass") && !treatmentApprochEntryObject.get("drugClass")
                                    .isJsonNull() ? extractDrugClass(treatmentApprochEntryObject.getAsJsonObject("drugClass")) : null)
                            .therapy(treatmentApprochEntryObject.has("therapy") && !treatmentApprochEntryObject.get("therapy").isJsonNull()
                                    ? extractTherapy(treatmentApprochEntryObject.getAsJsonObject("therapy"))
                                    : null)
                            .reference(extractReference(treatmentApprochEntryObject.getAsJsonArray("references")))
                            .createDate(DateConverter.convertDate(JsonFunctions.string(treatmentApprochEntryObject, "createDate")))
                            .updateDate(DateConverter.convertDate(JsonFunctions.string(treatmentApprochEntryObject, "updateDate")))
                            .build());
                }
                reader.close();
            }
        }
        LOGGER.info("Finished reading treatment approach dir");

        return treatmentAppraoch;
    }

    @NotNull
    private static DrugClassInfo extractDrugClass(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker drugClassObjectChecker = TreatmentApprochDataModelChecker.drugClassObjectChecker();
        drugClassObjectChecker.check(jsonObject);

        return ImmutableDrugClassInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .drugClass(JsonFunctions.string(jsonObject, "drugClass"))
                .build();
    }

    @NotNull
    private static TherapyInfo extractTherapy(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker therapyObjectChecker = TreatmentApprochDataModelChecker.therapyObjectChecker();
        therapyObjectChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.nullableString(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    private static List<ReferenceInfo> extractReference(@NotNull JsonArray jsonArray) {
        List<ReferenceInfo> references = Lists.newArrayList();
        JsonDatamodelChecker referenceObjectChecker = TreatmentApprochDataModelChecker.referenceObjectChecker();

        for (JsonElement reference : jsonArray) {
            JsonObject referenceJsonObject = reference.getAsJsonObject();
            referenceObjectChecker.check(referenceJsonObject);

            references.add(ImmutableReferenceInfo.builder()
                    .id(JsonFunctions.integer(referenceJsonObject, "id"))
                    .pubMedId(JsonFunctions.nullableString(referenceJsonObject, "pubMedId"))
                    .title(JsonFunctions.string(referenceJsonObject, "title"))
                    .url(JsonFunctions.nullableString(referenceJsonObject, "url"))
                    .build());
        }
        return references;
    }
}
