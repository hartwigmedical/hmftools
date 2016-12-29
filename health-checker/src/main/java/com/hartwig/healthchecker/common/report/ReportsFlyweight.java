package com.hartwig.healthchecker.common.report;

import java.util.HashMap;
import java.util.Map;

import org.jetbrains.annotations.NotNull;

final class ReportsFlyweight {

    private static final String STDOUT = "stdout";
    private static final String JSON = "json";

    private static final Map<String, Report> FLYWEIGHT = new HashMap<>();
    private static final ReportsFlyweight INSTANCE = new ReportsFlyweight();

    static {
        FLYWEIGHT.put(JSON, JsonReport.getInstance());
        FLYWEIGHT.put(STDOUT, StandardOutputReport.getInstance());
    }

    private ReportsFlyweight() {
    }

    @NotNull
    static ReportsFlyweight getInstance() {
        return INSTANCE;
    }

    @NotNull
    Report getReport(@NotNull final String reportType) {
        Report defaultReport = StandardOutputReport.getInstance();

        if (FLYWEIGHT.containsKey(reportType)) {
            defaultReport = FLYWEIGHT.get(reportType);
        }
        return defaultReport;
    }
}
