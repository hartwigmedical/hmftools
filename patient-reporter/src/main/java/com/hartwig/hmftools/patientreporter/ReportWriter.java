package com.hartwig.hmftools.patientreporter;

import java.io.IOException;

import org.jetbrains.annotations.NotNull;

public interface ReportWriter {

    void writeSequenceReport(@NotNull AnalysedPatientReport report, @NotNull String outputFilePath) throws IOException;

    void writeNonSequenceableReport(@NotNull NotAnalysedPatientReport report, @NotNull String outputFilePath) throws IOException;
}
