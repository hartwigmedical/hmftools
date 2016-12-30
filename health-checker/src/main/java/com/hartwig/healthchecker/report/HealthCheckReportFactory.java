package com.hartwig.healthchecker.report;

import org.jetbrains.annotations.NotNull;

public final class HealthCheckReportFactory {

    private HealthCheckReportFactory() {
    }

    @NotNull
    public static Report create(@NotNull String reportType) {
        return ReportsFlyweight.getInstance().getReport(reportType);
    }
}
