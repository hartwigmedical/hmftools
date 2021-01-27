package com.hartwig.hmftools.ckb.reader.drugclass;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;

import org.jetbrains.annotations.NotNull;

public final class DrugClassDataModelChecker {

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
    public static JsonDatamodelChecker drugClassDrugsObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("drugName", true);
        map.put("terms", true);

        return new JsonDatamodelChecker("DrugClassDrugsObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugClassTreatmentApproachesObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("profileName", true);

        return new JsonDatamodelChecker("DrugClassTreatmentApproachesObject", map);
    }
}
