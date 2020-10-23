package com.hartwig.hmftools.patientdb.diseaseontology;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class JsonFunctions {
    private static final Logger LOGGER = LogManager.getLogger(JsonFunctions.class);

    @NotNull
    static String toString(@NotNull JsonObject object, @NotNull String field) {
        assert object.has(field);

        JsonElement element = object.get(field);
        if (!element.isJsonPrimitive()) {
            LOGGER.warn("Converting {} to String for element {}.", field, element);
        }
        return element.getAsJsonPrimitive().getAsString();
    }

    @NotNull
    static List<String> stringList(@NotNull JsonObject object, @NotNull String field) {
        assert object.has(field);

        if (object.get(field).isJsonNull()) {
            return Lists.newArrayList();
        }

        List<String> values = Lists.newArrayList();
        if (object.get(field).isJsonPrimitive()) {
            values.add(toString(object, field));
        } else {
            assert object.get(field).isJsonArray();
            for (JsonElement element : object.getAsJsonArray(field)) {
                if (!element.isJsonPrimitive()) {
                    LOGGER.warn("Converting array value for {} into string for element {}", field, element);
                }
                values.add(element.getAsJsonPrimitive().getAsString());
            }
        }
        return values;
    }

    @NotNull
    static List<String> optionalStringList(@NotNull JsonObject object, @NotNull String field) {
        if (!object.has(field)) {
            return Lists.newArrayList();
        }

        return stringList(object, field);
    }

}
