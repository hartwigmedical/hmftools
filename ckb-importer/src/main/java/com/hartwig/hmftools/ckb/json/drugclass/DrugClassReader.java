package com.hartwig.hmftools.ckb.json.drugclass;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.ckb.json.CkbJsonDirectoryReader;
import com.hartwig.hmftools.ckb.json.common.DrugInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableDrugInfo;
import com.hartwig.hmftools.ckb.json.common.ImmutableTreatmentApproachInfo;
import com.hartwig.hmftools.ckb.json.common.TreatmentApproachInfo;
import com.hartwig.hmftools.ckb.util.DateConverter;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;
import com.hartwig.hmftools.common.utils.json.JsonFunctions;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class DrugClassReader extends CkbJsonDirectoryReader<DrugClass> {

    public DrugClassReader(@Nullable final Integer maxFilesToRead) {
        super(maxFilesToRead);
    }

    @NotNull
    @Override
    protected DrugClass read(@NotNull final JsonObject object) {
        JsonDatamodelChecker drugsClassChecker = DrugClassDataModelChecker.drugClassObjectChecker();
        drugsClassChecker.check(object);

        return ImmutableDrugClass.builder()
                .id(JsonFunctions.integer(object, "id"))
                .drugClass(JsonFunctions.string(object, "drugClass"))
                .createDate(DateConverter.toDate(JsonFunctions.string(object, "createDate")))
                .drugs(extractDrugs(object.getAsJsonArray("drugs")))
                .treatmentApproaches(extractTreatmentApproaches(object.getAsJsonArray("treatmentApproaches")))
                .build();
    }

    @NotNull
    private static List<DrugInfo> extractDrugs(@NotNull JsonArray jsonArray) {
        List<DrugInfo> drugs = Lists.newArrayList();
        JsonDatamodelChecker drugChecker = DrugClassDataModelChecker.drugObjectChecker();

        for (JsonElement drug : jsonArray) {
            JsonObject drugObject = drug.getAsJsonObject();
            drugChecker.check(drugObject);

            drugs.add(ImmutableDrugInfo.builder()
                    .id(JsonFunctions.integer(drugObject, "id"))
                    .drugName(JsonFunctions.string(drugObject, "drugName"))
                    .terms(JsonFunctions.optionalStringList(drugObject, "terms"))
                    .build());
        }
        return drugs;
    }

    @NotNull
    private static List<TreatmentApproachInfo> extractTreatmentApproaches(@NotNull JsonArray jsonArray) {
        List<TreatmentApproachInfo> treatmentApproaches = Lists.newArrayList();
        JsonDatamodelChecker treatmentApproachChecker = DrugClassDataModelChecker.treatmentApproachObjectChecker();

        for (JsonElement treatmentApproach : jsonArray) {
            JsonObject treatmentApproachObject = treatmentApproach.getAsJsonObject();
            treatmentApproachChecker.check(treatmentApproachObject);

            treatmentApproaches.add(ImmutableTreatmentApproachInfo.builder()
                    .id(JsonFunctions.integer(treatmentApproachObject, "id"))
                    .name(JsonFunctions.string(treatmentApproachObject, "name"))
                    .profileName(JsonFunctions.string(treatmentApproachObject, "profileName"))
                    .build());
        }
        return treatmentApproaches;
    }
}
