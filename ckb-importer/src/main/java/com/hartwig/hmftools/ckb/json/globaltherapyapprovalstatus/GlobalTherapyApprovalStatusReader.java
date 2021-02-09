package com.hartwig.hmftools.ckb.json.globaltherapyapprovalstatus;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.ckb.json.CkbJsonDirectoryReader;
import com.hartwig.hmftools.ckb.json.common.GlobalApprovalStatusInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableGlobalApprovalStatusInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableIndicationInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableMolecularProfileInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableTherapyInfo;
import com.hartwig.hmftools.ckb.json.common.IndicationInfo;
import com.hartwig.hmftools.ckb.json.common.MolecularProfileInfo;
import com.hartwig.hmftools.ckb.json.common.TherapyInfo;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class GlobalTherapyApprovalStatusReader extends CkbJsonDirectoryReader<GlobalTherapyApprovalStatus> {

    public GlobalTherapyApprovalStatusReader(@Nullable final Integer maxFilesToRead) {
        super(maxFilesToRead);
    }

    @NotNull
    @Override
    protected GlobalTherapyApprovalStatus read(@NotNull final JsonObject object) {
        JsonDatamodelChecker globalTherapyApprovalStatusChecker =
                GlobalTherapyApprovalStatusDataModelChecker.globalTherapyApprovalStatusObjectChecker();
        globalTherapyApprovalStatusChecker.check(object);

        return ImmutableGlobalTherapyApprovalStatus.builder()
                .totalCount(JsonFunctions.integer(object, "totalCount"))
                .globalApprovalStatus(extractGlobalTherapyApprovalStatus(object.getAsJsonArray("globalTherapyApprovalStatuses")))
                .build();
    }

    @NotNull
    private static List<GlobalApprovalStatusInfo> extractGlobalTherapyApprovalStatus(@NotNull JsonArray jsonArray) {
        List<GlobalApprovalStatusInfo> globalTherapyApprovalStatusses = Lists.newArrayList();
        JsonDatamodelChecker globalTherapyApprovalStatusListChecker =
                GlobalTherapyApprovalStatusDataModelChecker.globalTherapyApprovalStatusListObjectChecker();

        for (JsonElement globalTherapyApprovalStatus : jsonArray) {
            JsonObject globalTherapyApprovalStatusJsonObject = globalTherapyApprovalStatus.getAsJsonObject();
            globalTherapyApprovalStatusListChecker.check(globalTherapyApprovalStatusJsonObject);

            globalTherapyApprovalStatusses.add(ImmutableGlobalApprovalStatusInfo.builder()
                    .id(JsonFunctions.integer(globalTherapyApprovalStatusJsonObject, "id"))
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
    private static TherapyInfo extractGlobalTherapyApprovalStatusTherapy(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker globalTherapyApprovalStatusTherapyChecker =
                GlobalTherapyApprovalStatusDataModelChecker.globalTherapyApprovalStatusTherapyObjectChecker();
        globalTherapyApprovalStatusTherapyChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.nullableString(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    private static IndicationInfo extractGlobalTherapyApprovalStatusIndication(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker globalTherapyApprovalStatusIndicationChecker =
                GlobalTherapyApprovalStatusDataModelChecker.globalTherapyApprovalStatusIndicationObjectChecker();
        globalTherapyApprovalStatusIndicationChecker.check(jsonObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    private static MolecularProfileInfo extractGlobalTherapyApprovalStatusMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker globalTherapyApprovalStatusMolecularProfileChecker =
                GlobalTherapyApprovalStatusDataModelChecker.globalTherapyApprovalStatusMolecularProfileObjectChecker();
        globalTherapyApprovalStatusMolecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }
}
