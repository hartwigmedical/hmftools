package com.hartwig.hmftools.patientreporter.report;

import java.io.IOException;

import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.PatientReport;
import com.hartwig.hmftools.patientreporter.algo.NotSequenceableReason;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.exception.DRException;

public interface ReportWriter {

    void writeSequenceReport(@NotNull final PatientReport report, @NotNull final HmfReporterData reporterData)
            throws IOException, DRException;

    void writeNonSequenceableReport(@NotNull final String sample, @NotNull final String tumorType,
            @NotNull final String tumorPercentage, @NotNull final NotSequenceableReason reason)
            throws IOException, DRException;
}
