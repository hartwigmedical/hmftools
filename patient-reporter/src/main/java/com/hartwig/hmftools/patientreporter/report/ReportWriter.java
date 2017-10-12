package com.hartwig.hmftools.patientreporter.report;

import java.io.IOException;

import com.hartwig.hmftools.patientreporter.HmfReporterData;
import com.hartwig.hmftools.patientreporter.NotSequencedPatientReport;
import com.hartwig.hmftools.patientreporter.SequencedPatientReport;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.exception.DRException;

public interface ReportWriter {

    void writeSequenceReport(@NotNull final SequencedPatientReport report, @NotNull final HmfReporterData reporterData)
            throws IOException, DRException;

    void writeNonSequenceableReport(@NotNull final NotSequencedPatientReport report) throws IOException, DRException;
}
