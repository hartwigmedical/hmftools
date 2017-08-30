package com.hartwig.hmftools.healthchecker.report;

import java.util.Optional;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.common.context.RunContext;
import com.hartwig.hmftools.common.exception.GenerateReportException;

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
    public Optional<String> generateReport(@NotNull final RunContext runContext, @NotNull final String fileName)
            throws GenerateReportException {
        final JsonArray reportArray = computeElements();

        final JsonObject reportJson = new JsonObject();
        reportJson.add("health_checks", reportArray);

        return Optional.ofNullable(AbstractJsonBaseReport.GSON.toJson(reportJson));
    }
}
