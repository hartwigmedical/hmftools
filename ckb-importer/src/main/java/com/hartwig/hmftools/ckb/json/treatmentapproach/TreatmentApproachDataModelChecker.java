package com.hartwig.hmftools.ckb.json.treatmentapproach;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;

import org.jetbrains.annotations.NotNull;

final class TreatmentApproachDataModelChecker {

    private TreatmentApproachDataModelChecker() {
    }

    @NotNull
    public static JsonDatamodelChecker treatmentApproachObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("profileName", true);
        map.put("drugClass", true);
        map.put("therapy", true);
        map.put("references", true);
        map.put("createDate", true);
        map.put("updateDate", true);

        return new JsonDatamodelChecker("TreatmentApproachObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugClassObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("drugClass", true);

        return new JsonDatamodelChecker("TreatmentApproachDrugClassObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker therapyObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapyName", true);
        map.put("synonyms", true);

        return new JsonDatamodelChecker("TreatmentApproachTherapyObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker referenceObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("pubMedId", true);
        map.put("title", true);
        map.put("url", true);

        return new JsonDatamodelChecker("TreatmentApproachReferenceObject", map);
    }
}
