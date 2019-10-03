package com.hartwig.hmftools.vicc.reader;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;

import org.jetbrains.annotations.NotNull;

final class JsonFunctions {

    private JsonFunctions() {
    }

    @NotNull
    static List<String> jsonArrayToStringList(@NotNull JsonArray array) {
        List<String> values = Lists.newArrayList();
        for (JsonElement element : array) {
            values.add(element.getAsString());
        }
        return values;
    }
}
