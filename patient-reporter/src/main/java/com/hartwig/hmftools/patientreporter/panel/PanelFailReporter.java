package com.hartwig.hmftools.patientreporter.panel;

import com.hartwig.hmftools.patientreporter.qcfail.QCFailReporter;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class PanelFailReporter {

    private static final Logger LOGGER = LogManager.getLogger(QCFailReporter.class);

    @NotNull
    private final QCFailPanelReportData reportData;
    private final String reportDate;

    public PanelFailReporter(@NotNull final QCFailPanelReportData reportData, @NotNull final String reportDate) {
        this.reportData = reportData;
        this.reportDate = reportDate;
    }

    @NotNull
    public PanelFailReport run() {
        return ImmutablePanelFailReport.builder().build();
    }
}