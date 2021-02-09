package com.hartwig.hmftools.ckb.reader.treatmentapproch;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.ckb.datamodel.common.DrugClassInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableDrugClassInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableReferenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ImmutableTherapyInfo;
import com.hartwig.hmftools.ckb.datamodel.common.ReferenceInfo;
import com.hartwig.hmftools.ckb.datamodel.common.TherapyInfo;
import com.hartwig.hmftools.ckb.datamodel.treatmentapproach.ImmutableTreatmentApproach;
import com.hartwig.hmftools.ckb.datamodel.treatmentapproach.TreatmentApproach;
import com.hartwig.hmftools.ckb.reader.CkbJsonDirectoryReader;
import com.hartwig.hmftools.ckb.util.DateConverter;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class TreatmentApproachReader extends CkbJsonDirectoryReader<TreatmentApproach> {

    public TreatmentApproachReader(@Nullable final Integer maxFilesToRead) {
        super(maxFilesToRead);
    }

    @NotNull
    @Override
    protected TreatmentApproach read(@NotNull final JsonObject object) {
        JsonDatamodelChecker treatmentApprochObjectChecker = TreatmentApproachDataModelChecker.treatmentApprochObjectChecker();
        treatmentApprochObjectChecker.check(object);

        return ImmutableTreatmentApproach.builder()
                .id(JsonFunctions.integer(object, "id"))
                .name(JsonFunctions.string(object, "name"))
                .profileName(JsonFunctions.string(object, "profileName"))
                .drugClass(object.has("drugClass") && !object.get("drugClass").isJsonNull() ? extractDrugClass(object.getAsJsonObject(
                        "drugClass")) : null)
                .therapy(object.has("therapy") && !object.get("therapy").isJsonNull()
                        ? extractTherapy(object.getAsJsonObject("therapy"))
                        : null)
                .reference(extractReference(object.getAsJsonArray("references")))
                .createDate(DateConverter.toDate(JsonFunctions.string(object, "createDate")))
                .updateDate(DateConverter.toDate(JsonFunctions.string(object, "updateDate")))
                .build();
    }

    @NotNull
    private static DrugClassInfo extractDrugClass(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker drugClassObjectChecker = TreatmentApproachDataModelChecker.drugClassObjectChecker();
        drugClassObjectChecker.check(jsonObject);

        return ImmutableDrugClassInfo.builder()
                .id(JsonFunctions.integer(jsonObject, "id"))
                .drugClass(JsonFunctions.string(jsonObject, "drugClass"))
                .build();
    }

    @NotNull
    private static TherapyInfo extractTherapy(@NotNull JsonObject jsonObject) {
        JsonDatamodelChecker therapyObjectChecker = TreatmentApproachDataModelChecker.therapyObjectChecker();
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
        JsonDatamodelChecker referenceObjectChecker = TreatmentApproachDataModelChecker.referenceObjectChecker();

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
