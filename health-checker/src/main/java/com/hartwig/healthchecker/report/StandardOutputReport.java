package com.hartwig.healthchecker.report;

import java.util.Optional;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.hartwig.healthchecker.exception.GenerateReportException;
import com.hartwig.healthchecker.io.dir.RunContext;

import org.jetbrains.annotations.NotNull;

final class StandardOutputReport extends AbstractJsonBaseReport {

    private static final StandardOutputReport INSTANCE = new StandardOutputReport();

    private StandardOutputReport() {
    }

    static StandardOutputReport getInstance() {
        return INSTANCE;
    }

    @NotNull
    @Override
    public Optional<String> generateReport(@NotNull final RunContext runContext, @NotNull final String outputPath)
            throws GenerateReportException {
        final JsonArray reportArray = computeElements();

        final JsonObject reportJson = new JsonObject();
        reportJson.add("health_checks", reportArray);

        return Optional.ofNullable(AbstractJsonBaseReport.GSON.toJson(reportJson));
    }
}
