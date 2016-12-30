package com.hartwig.healthchecker.report;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

import com.google.gson.FieldNamingPolicy;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.healthchecker.runners.CheckType;
import com.hartwig.healthchecker.result.BaseResult;

import org.jetbrains.annotations.NotNull;

abstract class AbstractJsonBaseReport implements Report {

    static final Gson GSON = new GsonBuilder().setPrettyPrinting()
            .setFieldNamingPolicy(FieldNamingPolicy.LOWER_CASE_WITH_UNDERSCORES).disableHtmlEscaping().create();

    private static final Map<CheckType, BaseResult> HEALTH_CHECKS = new ConcurrentHashMap<>();

    @Override
    public void addResult(@NotNull final BaseResult result) {
        HEALTH_CHECKS.putIfAbsent(result.getCheckType(), result);
    }

    @NotNull
    JsonArray computeElements() {
        final JsonArray reportArray = new JsonArray();

        HEALTH_CHECKS.forEach((checkType, baseReport) -> {
            final JsonElement configJson = GSON.toJsonTree(baseReport);

            final JsonObject element = new JsonObject();
            element.add(checkType.toString(), configJson);

            reportArray.add(element);
        });

        return reportArray;
    }
}
