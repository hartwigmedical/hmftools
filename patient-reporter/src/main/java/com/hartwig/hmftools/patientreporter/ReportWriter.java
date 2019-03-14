package com.hartwig.hmftools.patientreporter;

import java.io.IOException;

import org.jetbrains.annotations.NotNull;

public interface ReportWriter {

    void writeAnalysedPatientReport(@NotNull AnalysedPatientReport report, @NotNull String outputFilePath) throws IOException;

    void writeQCFailReport(@NotNull QCFailReport report, @NotNull String outputFilePath) throws IOException;
}
