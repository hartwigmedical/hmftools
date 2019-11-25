package com.hartwig.hmftools.vicc.reader;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

final class JsonFunctions {

    private JsonFunctions() {
    }

    @NotNull
    static List<String> optionalStringList(@NotNull JsonObject object, @NotNull String field) {
        if (!object.has(field)) {
            return Lists.newArrayList();
        }

        assert object.get(field).isJsonArray();
        return jsonArrayToStringList(object.getAsJsonArray(field));
    }

    @NotNull
    static List<String> jsonArrayToStringList(@NotNull JsonArray array) {
        List<String> values = Lists.newArrayList();
        for (JsonElement element : array) {
            values.add(element.getAsString());
        }
        return values;
    }

    @Nullable
    static String optionalString(@NotNull JsonObject object, @NotNull String field) {
        if (!object.has(field)) {
            return null;
        }

        JsonElement element = object.get(field);
        assert element.isJsonPrimitive();
        return element.getAsJsonPrimitive().getAsString();
    }

    @Nullable
    static String nullableString(@NotNull JsonObject object, @NotNull String field) {
        assert object.has(field);
        JsonElement element = object.get(field);

        if (element.isJsonNull()) {
            return null;
        }

        assert element.isJsonPrimitive();
        return element.getAsJsonPrimitive().getAsString();
    }

    @Nullable
    static String optionalNullableString(@NotNull JsonObject object, @NotNull String field) {
        if (!object.has(field)) {
            return null;
        }

        JsonElement element = object.get(field);

        if (element.isJsonNull()) {
            return null;
        }

        assert element.isJsonPrimitive();
        return element.getAsJsonPrimitive().getAsString();
    }
}
