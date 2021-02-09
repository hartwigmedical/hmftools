package com.hartwig.hmftools.ckb.json.drugclass;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;

import org.jetbrains.annotations.NotNull;

final class DrugClassDataModelChecker {

    private DrugClassDataModelChecker(){
    }

    @NotNull
    public static JsonDatamodelChecker drugClassObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("drugClass", true);
        map.put("createDate", true);
        map.put("drugs", true);
        map.put("treatmentApproaches", true);

        return new JsonDatamodelChecker("DrugClassObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("drugName", true);
        map.put("terms", true);

        return new JsonDatamodelChecker("DrugClassDrugObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker treatmentApproachObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("profileName", true);

        return new JsonDatamodelChecker("DrugClassTreatmentApproachObject", map);
    }
}
