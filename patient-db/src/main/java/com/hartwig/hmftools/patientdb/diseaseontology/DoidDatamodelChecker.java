package com.hartwig.hmftools.patientdb.diseaseontology;

import java.util.Map;
import java.util.Set;

import com.google.gson.JsonObject;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class DoidDatamodelChecker {

    private static final Logger LOGGER = LogManager.getLogger(DoidDatamodelChecker.class);

    @NotNull
    private final String name;
    @NotNull
    private final Map<String, Boolean> datamodel;

    DoidDatamodelChecker(@NotNull final String name, @NotNull final Map<String, Boolean> datamodel) {
        this.name = name;
        this.datamodel = datamodel;
    }

    void check(@NotNull JsonObject object) {
        Set<String> keys = object.keySet();

        for (String key : keys) {
            if (!datamodel.containsKey(key)) {
                LOGGER.warn("JSON object contains key '{}' which is not expected for '{}'", key, name);
            }
        }

        for (Map.Entry<String, Boolean> datamodelEntry : datamodel.entrySet()) {
            if (datamodelEntry.getValue() && !keys.contains(datamodelEntry.getKey())) {
                LOGGER.warn("Mandatory key '{}' missing from JSON object in '{}'", datamodelEntry.getKey(), name);
            }
        }
    }
}
