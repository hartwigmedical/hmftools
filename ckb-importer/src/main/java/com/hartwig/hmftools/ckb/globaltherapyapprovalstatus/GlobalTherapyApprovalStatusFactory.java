package com.hartwig.hmftools.ckb.globaltherapyapprovalstatus;

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

public class GlobalTherapyApprovalStatusFactory {

    private static final Logger LOGGER = LogManager.getLogger(GlobalTherapyApprovalStatusFactory.class);

    private GlobalTherapyApprovalStatusFactory() {

    }

    @NotNull
    public static List<GlobalTherapyApprovalStatus> readingGlobalTherapyApprovalStatus(@NotNull String globalTherapyApprovalStatusDir)
            throws IOException {
        LOGGER.info("Start reading genes");

        List<GlobalTherapyApprovalStatus> globalTherapyApprovalStatusses = Lists.newArrayList();
        File[] filesGlobalTherapyApprovalStatus = new File(globalTherapyApprovalStatusDir).listFiles();

        if (filesGlobalTherapyApprovalStatus != null) {
            LOGGER.info("The total files in the genes dir is {}", filesGlobalTherapyApprovalStatus.length);

            for (File globalTherapyApprovalStatus : filesGlobalTherapyApprovalStatus) {
                JsonParser parser = new JsonParser();
                JsonReader reader = new JsonReader(new FileReader(globalTherapyApprovalStatus));
                reader.setLenient(true);

                while (reader.peek() != JsonToken.END_DOCUMENT) {
                    JsonObject globalTherapyApprovalStatusEntryObject = parser.parse(reader).getAsJsonObject();
                    JsonDatamodelChecker globalTherapyApprovalStatusChecker =
                            GlobalTherapyApprovalStatusDataModelChecker.globalTherapyApprovalStatusObjectChecker();
                    globalTherapyApprovalStatusChecker.check(globalTherapyApprovalStatusEntryObject);

                    globalTherapyApprovalStatusses.add(ImmutableGlobalTherapyApprovalStatus.builder()
                            .totalCount(JsonFunctions.string(globalTherapyApprovalStatusEntryObject, "totalCount"))
                            .approvalStatusList(extractGlobalTherapyApprovalStatus(globalTherapyApprovalStatusEntryObject.getAsJsonArray(
                                    "globalTherapyApprovalStatuses")))
                            .build());
                }
            }
        }
        return globalTherapyApprovalStatusses;
    }

    @NotNull
    public static List<GlobalTherapyApprovalStatusList> extractGlobalTherapyApprovalStatus(@NotNull JsonArray jsonArray) {
        List<GlobalTherapyApprovalStatusList> globalTherapyApprovalStatusses = Lists.newArrayList();
        JsonDatamodelChecker globalTherapyApprovalStatusListChecker =
                GlobalTherapyApprovalStatusDataModelChecker.globalTherapyApprovalStatusListObjectChecker();

        for (JsonElement globalTherapyApprovalStatus : jsonArray) {
            JsonObject globalTherapyApprovalStatusJsonObject = globalTherapyApprovalStatus.getAsJsonObject();
            globalTherapyApprovalStatusListChecker.check(globalTherapyApprovalStatusJsonObject);

            globalTherapyApprovalStatusses.add(ImmutableGlobalTherapyApprovalStatusList.builder()
                    .id(JsonFunctions.string(globalTherapyApprovalStatusJsonObject, "id"))
                    .therapy(extractGlobalTherapyApprovalStatusTherapy(globalTherapyApprovalStatusJsonObject.getAsJsonObject("therapy")))
                    .indication(extractGlobalTherapyApprovalStatusIndication(globalTherapyApprovalStatusJsonObject.getAsJsonObject(
                            "indication")))
                    .molecularProfile(extractGlobalTherapyApprovalStatusMolecularProfile(globalTherapyApprovalStatusJsonObject.getAsJsonObject(
                            "molecularProfile")))
                    .approvalAuthority(JsonFunctions.string(globalTherapyApprovalStatusJsonObject, "approvalAuthority"))
                    .approvalStatus(JsonFunctions.string(globalTherapyApprovalStatusJsonObject, "approvalStatus"))
                    .build());
        }
        return globalTherapyApprovalStatusses;
    }

    @NotNull
    public static GlobalTherapyApprovalStatusTherapy extractGlobalTherapyApprovalStatusTherapy(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker globalTherapyApprovalStatusTherapyChecker =
                GlobalTherapyApprovalStatusDataModelChecker.globalTherapyApprovalStatusTherapyObjectChecker();
        globalTherapyApprovalStatusTherapyChecker.check(jsonObject);

        return ImmutableGlobalTherapyApprovalStatusTherapy.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.nullableString(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    public static GlobalTherapyApprovalStatusIndication extractGlobalTherapyApprovalStatusIndication(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker globalTherapyApprovalStatusIndicationChecker =
                GlobalTherapyApprovalStatusDataModelChecker.globalTherapyApprovalStatusIndicationObjectChecker();
        globalTherapyApprovalStatusIndicationChecker.check(jsonObject);

        return ImmutableGlobalTherapyApprovalStatusIndication.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    public static GlobalTherapyApprovalStatusMolecularProfile extractGlobalTherapyApprovalStatusMolecularProfile(
            @NotNull JsonObject jsonObject) {
        JsonDatamodelChecker globalTherapyApprovalStatusMolecularProfileChecker =
                GlobalTherapyApprovalStatusDataModelChecker.globalTherapyApprovalStatusMolecularprofileObjectChecker();
        globalTherapyApprovalStatusMolecularProfileChecker.check(jsonObject);

        return ImmutableGlobalTherapyApprovalStatusMolecularProfile.builder()
                .id(JsonFunctions.string(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }

}
