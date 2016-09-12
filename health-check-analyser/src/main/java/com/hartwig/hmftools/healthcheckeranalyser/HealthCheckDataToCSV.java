package com.hartwig.hmftools.healthcheckeranalyser;

import com.hartwig.hmftools.healthcheckeranalyser.model.HealthCheckReport;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class HealthCheckDataToCSV {

    private HealthCheckDataToCSV() {
    }

    @NotNull
    public static String header(@NotNull HealthCheckReport report) {
        return Strings.EMPTY;
    }
}
