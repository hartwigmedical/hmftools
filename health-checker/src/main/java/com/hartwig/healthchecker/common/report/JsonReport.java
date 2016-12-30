package com.hartwig.healthchecker.common.report;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Optional;

import com.google.gson.JsonArray;
import com.google.gson.JsonObject;
import com.hartwig.healthchecker.common.exception.GenerateReportException;
import com.hartwig.healthchecker.common.io.dir.RunContext;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class JsonReport extends AbstractJsonBaseReport {

    private static final JsonReport INSTANCE = new JsonReport();

    private static final Logger LOGGER = LogManager.getLogger(JsonReport.class);
    private static final String REPORT_NAME = "%s_health_checks_%s.json";

    private static final String ERROR_GENERATING_REPORT = "Error occurred whilst generating reports. Error -> %s";

    private JsonReport() {
    }

    static JsonReport getInstance() {
        return INSTANCE;
    }

    @NotNull
    @Override
    public Optional<String> generateReport(@NotNull final RunContext runContext, @NotNull final String outputPath)
            throws GenerateReportException {
        final JsonArray reportArray = computeElements();

        final String runName = toName(runContext.runDirectory());
        final String fileName = String.format("%s/%s", outputPath,
                String.format(REPORT_NAME, runName, System.currentTimeMillis()));

        try (FileWriter fileWriter = new FileWriter(new File(fileName))) {
            final JsonObject reportJson = new JsonObject();
            reportJson.add("health_checks", reportArray);
            fileWriter.write(AbstractJsonBaseReport.GSON.toJson(reportJson));
            fileWriter.flush();
        } catch (final IOException e) {
            LOGGER.error(String.format(ERROR_GENERATING_REPORT, e.getMessage()));
            throw (GenerateReportException) new GenerateReportException(e.getMessage()).initCause(e);
        }
        return Optional.of(fileName);
    }

    @NotNull
    private static String toName(@NotNull final String runDirectory) {
        final String[] parts = runDirectory.split(File.separator);
        return parts[parts.length - 1];
    }
}
