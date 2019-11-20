package com.hartwig.hmftools.vicc.reader;

import java.util.Map;
import java.util.Set;

import com.google.gson.JsonObject;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

class ViccDatamodelChecker {

    private static final Logger LOGGER = LogManager.getLogger(ViccDatamodelChecker.class);

    @NotNull
    private final String name;
    @NotNull
    private final Map<String, Boolean> datamodel;

    ViccDatamodelChecker(@NotNull final String name, @NotNull final Map<String, Boolean> datamodel) {
        this.name = name;
        this.datamodel = datamodel;
    }

    boolean check(@NotNull JsonObject object) {
        boolean correct = true;
        Set<String> keys = object.keySet();

        for (String key : keys) {
            if (!datamodel.containsKey(key)) {
                LOGGER.warn("JSON object contains key '{}' which is not expected for '{}'", key, name);
                correct = false;
            }
        }

        for (Map.Entry<String, Boolean> datamodelEntry : datamodel.entrySet()) {
            Boolean isMandatory = datamodelEntry.getValue();
            assert isMandatory != null;
            if (isMandatory && !keys.contains(datamodelEntry.getKey())) {
                LOGGER.warn("Mandatory key '{}' missing from JSON object in '{}'", datamodelEntry.getKey(), name);
                correct = false;
            }
        }

        return correct;
    }
}
