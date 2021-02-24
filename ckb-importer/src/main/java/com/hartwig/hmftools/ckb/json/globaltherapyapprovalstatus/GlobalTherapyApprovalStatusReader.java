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

public class GlobalTherapyApprovalStatusReader extends CkbJsonDirectoryReader<JsonGlobalTherapyApprovalStatus> {

    public GlobalTherapyApprovalStatusReader(@Nullable final Integer maxFilesToRead) {
        super(maxFilesToRead);
    }

    @NotNull
    @Override
    protected JsonGlobalTherapyApprovalStatus read(@NotNull final JsonObject object) {
        JsonDatamodelChecker statusChecker = GlobalTherapyApprovalStatusDataModelChecker.globalTherapyApprovalStatusObjectChecker();
        statusChecker.check(object);

        return ImmutableJsonGlobalTherapyApprovalStatus.builder()
                .totalCount(JsonFunctions.integer(object, "totalCount"))
                .globalApprovalStatuses(extractStatuses(object.getAsJsonArray("globalTherapyApprovalStatuses")))
                .build();
    }

    @NotNull
    private static List<GlobalApprovalStatusInfo> extractStatuses(@NotNull JsonArray jsonArray) {
        List<GlobalApprovalStatusInfo> statuses = Lists.newArrayList();
        JsonDatamodelChecker listChecker = GlobalTherapyApprovalStatusDataModelChecker.listObjectChecker();

        for (JsonElement status : jsonArray) {
            JsonObject statusJsonObject = status.getAsJsonObject();
            listChecker.check(statusJsonObject);

            statuses.add(ImmutableGlobalApprovalStatusInfo.builder()
                    .id(JsonFunctions.integer(statusJsonObject, "id"))
                    .therapy(extractTherapy(statusJsonObject.getAsJsonObject("therapy")))
                    .indication(extractIndication(statusJsonObject.getAsJsonObject("indication")))
                    .molecularProfile(extractMolecularProfile(statusJsonObject.getAsJsonObject("molecularProfile")))
                    .approvalAuthority(JsonFunctions.string(statusJsonObject, "approvalAuthority"))
                    .approvalStatus(JsonFunctions.string(statusJsonObject, "approvalStatus"))
                    .build());
        }
        return statuses;
    }

    @NotNull
    private static TherapyInfo extractTherapy(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker therapyChecker = GlobalTherapyApprovalStatusDataModelChecker.therapyObjectChecker();
        therapyChecker.check(jsonObject);

        return ImmutableTherapyInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .therapyName(JsonFunctions.string(jsonObject, "therapyName"))
                .synonyms(JsonFunctions.optionalStringList(jsonObject, "synonyms"))
                .build();
    }

    @NotNull
    private static IndicationInfo extractIndication(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker indicationChecker = GlobalTherapyApprovalStatusDataModelChecker.indicationObjectChecker();
        indicationChecker.check(jsonObject);

        return ImmutableIndicationInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .name(JsonFunctions.string(jsonObject, "name"))
                .source(JsonFunctions.string(jsonObject, "source"))
                .build();
    }

    @NotNull
    private static MolecularProfileInfo extractMolecularProfile(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker molecularProfileChecker = GlobalTherapyApprovalStatusDataModelChecker.molecularProfileObjectChecker();
        molecularProfileChecker.check(jsonObject);

        return ImmutableMolecularProfileInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .profileName(JsonFunctions.string(jsonObject, "profileName"))
                .build();
    }
}
