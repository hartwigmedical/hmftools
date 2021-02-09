package com.hartwig.hmftools.ckb.json.treatmentapproch;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;

import org.jetbrains.annotations.NotNull;

final class TreatmentApproachDataModelChecker {

    private TreatmentApproachDataModelChecker() {
    }

    @NotNull
    public static JsonDatamodelChecker treatmentApprochObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("name", true);
        map.put("profileName", true);
        map.put("drugClass", true);
        map.put("therapy", true);
        map.put("references", true);
        map.put("createDate", true);
        map.put("updateDate", true);

        return new JsonDatamodelChecker("TreatmentApprochObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker drugClassObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("drugClass", true);

        return new JsonDatamodelChecker("TreatmentApprochDrugClassObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker therapyObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("therapyName", true);
        map.put("synonyms", true);

        return new JsonDatamodelChecker("TreatmentApprochTherapyObject", map);
    }

    @NotNull
    public static JsonDatamodelChecker referenceObjectChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("id", true);
        map.put("pubMedId", true);
        map.put("title", true);
        map.put("url", true);

        return new JsonDatamodelChecker("TreatmentApprochReferenceObject", map);
    }
}
