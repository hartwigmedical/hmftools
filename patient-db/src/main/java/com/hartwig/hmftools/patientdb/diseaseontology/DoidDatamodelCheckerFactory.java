package com.hartwig.hmftools.patientdb.diseaseontology;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.json.JsonDatamodelChecker;

import org.jetbrains.annotations.NotNull;

final class DoidDatamodelCheckerFactory {

    @NotNull
    static JsonDatamodelChecker doidBasicPropertyValuesChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("pred", true);
        map.put("val", true);
        return new JsonDatamodelChecker("doidBasicPropertyValues", map);
    }

    @NotNull
    static JsonDatamodelChecker doidMetadataChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("xrefs", false);
        map.put("synonyms", false);
        map.put("basicPropertyValues", false);
        map.put("definition", false);
        map.put("subsets", false);

        return new JsonDatamodelChecker("DoidMetadata", map);
    }
}
