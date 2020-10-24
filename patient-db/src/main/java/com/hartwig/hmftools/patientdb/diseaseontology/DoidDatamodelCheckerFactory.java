package com.hartwig.hmftools.patientdb.diseaseontology;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.utils.json.DatamodelChecker;

import org.jetbrains.annotations.NotNull;

final class DoidDatamodelCheckerFactory {

    @NotNull
    static DatamodelChecker doidBasicPropertyValuesChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("pred", true);
        map.put("val", true);
        return new DatamodelChecker("doidBasicPropertyValues", map);
    }

    @NotNull
    static DatamodelChecker doidMetadataChecker() {
        Map<String, Boolean> map = Maps.newHashMap();
        map.put("xrefs", false);
        map.put("synonyms", false);
        map.put("basicPropertyValues", false);
        map.put("definition", false);
        map.put("subsets", false);

        return new DatamodelChecker("DoidMetadata", map);
    }
}
