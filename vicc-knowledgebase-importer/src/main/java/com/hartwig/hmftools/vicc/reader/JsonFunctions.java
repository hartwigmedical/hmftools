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

    @Nullable
    static JsonObject optionalJsonObject(@NotNull JsonObject object, @NotNull String field) {
        if (!object.has(field)) {
            return null;
        }

        if (object.get(field).isJsonNull()) {
            return null;
        }

        assert object.get(field).isJsonObject();
        return object.getAsJsonObject(field);
    }

    @Nullable
    static JsonArray optionalJsonArray(@NotNull JsonObject object, @NotNull String field) {
        if (!object.has(field)) {
            return null;
        }

        if (object.get(field).isJsonNull()) {
            return null;
        }

        assert object.get(field).isJsonArray();
        return object.getAsJsonArray(field);
    }

    @NotNull
    static List<String> stringList(@NotNull JsonObject object, @NotNull String field) {
        assert object.has(field);

        if (object.get(field).isJsonNull()) {
            return Lists.newArrayList();
        }

        List<String> values = Lists.newArrayList();
        if (object.get(field).isJsonPrimitive()) {
            values.add(string(object, field));
        } else {
            assert object.get(field).isJsonArray();
            for (JsonElement element : object.getAsJsonArray(field)) {
                values.add(element.getAsString());
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

    @NotNull
    @Deprecated
    static List<String> toStringList(@NotNull JsonArray array) {
        List<String> values = Lists.newArrayList();
        for (JsonElement element : array) {
            values.add(element.getAsString());
        }
        return values;
    }

    @NotNull
    static String string(@NotNull JsonObject object, @NotNull String field) {
        assert object.has(field);

        JsonElement element = object.get(field);
        assert element.isJsonPrimitive();
        return element.getAsJsonPrimitive().getAsString();
    }

    @Nullable
    static String optionalString(@NotNull JsonObject object, @NotNull String field) {
        if (!object.has(field)) {
            return null;
        }

        return string(object, field);
    }

    @Nullable
    static String nullableString(@NotNull JsonObject object, @NotNull String field) {
        assert object.has(field);

        if (object.get(field).isJsonNull()) {
            return null;
        }

        return string(object, field);
    }

    @Nullable
    static String optionalNullableString(@NotNull JsonObject object, @NotNull String field) {
        if (!object.has(field)) {
            return null;
        }

        return nullableString(object, field);
    }
}
