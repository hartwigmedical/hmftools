package com.hartwig.hmftools.healthcheckeranalyser;

import java.io.FileReader;
import java.io.IOException;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.healthcheckeranalyser.model.HealthCheckReport;
import com.hartwig.hmftools.healthcheckeranalyser.model.HealthCheckReportFactory;

import org.jetbrains.annotations.NotNull;

final class HealthCheckReader {

    private static final Gson GSON = new GsonBuilder().create();

    private HealthCheckReader() {
    }

    @NotNull
    static HealthCheckReport readHealthCheckOutput(@NotNull String path) throws IOException {
        JsonObject json = GSON.fromJson(new FileReader(path), JsonObject.class);

        return HealthCheckReportFactory.fromHealthCheckReport(json);
    }
}
