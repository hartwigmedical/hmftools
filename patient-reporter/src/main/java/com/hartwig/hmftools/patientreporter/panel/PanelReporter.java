package com.hartwig.hmftools.patientreporter.panel;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PanelReporter {

    private static final Logger LOGGER = LogManager.getLogger(PanelReporter.class);

    @NotNull
    private final QCFailPanelReportData reportData;
    private final String reportDate;

    public PanelReporter(@NotNull final QCFailPanelReportData reportData, @NotNull final String reportDate) {
        this.reportData = reportData;
        this.reportDate = reportDate;
    }

    @NotNull
    public PanelReport run() {
        return ImmutablePanelReport.builder().build();
    }
}
