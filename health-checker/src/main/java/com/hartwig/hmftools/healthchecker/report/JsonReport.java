package com.hartwig.hmftools.healthchecker.report;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.Optional;

import com.google.common.collect.Maps;
import com.google.gson.FieldNamingPolicy;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.healthchecker.result.BaseResult;
import com.hartwig.hmftools.healthchecker.runners.CheckType;

import org.jetbrains.annotations.NotNull;

public class JsonReport implements Report {

    private static final Gson GSON = new GsonBuilder().setPrettyPrinting()
            .setFieldNamingPolicy(FieldNamingPolicy.LOWER_CASE_WITH_UNDERSCORES)
            .disableHtmlEscaping()
            .create();

    private final Map<CheckType, BaseResult> healthChecks = Maps.newHashMap();

    public JsonReport() {
    }

    @Override
    public void addResult(@NotNull final BaseResult result) {
        healthChecks.putIfAbsent(result.getCheckType(), result);
    }

    @NotNull
    @Override
    public Optional<String> generateReport(@NotNull final String fileName) throws IOException {
        final JsonArray reportArray = computeElements();

        FileWriter fileWriter = new FileWriter(new File(fileName));
        final JsonObject reportJson = new JsonObject();
        reportJson.add("health_checks", reportArray);

        fileWriter.write(GSON.toJson(reportJson));
        fileWriter.flush();

        return Optional.of(fileName);
    }

    @NotNull
    private JsonArray computeElements() {
        final JsonArray reportArray = new JsonArray();

        healthChecks.forEach((checkType, baseReport) -> {
            final JsonElement configJson = GSON.toJsonTree(baseReport);

            final JsonObject element = new JsonObject();
            element.add(checkType.toString(), configJson);

            reportArray.add(element);
        });

        return reportArray;
    }
}
