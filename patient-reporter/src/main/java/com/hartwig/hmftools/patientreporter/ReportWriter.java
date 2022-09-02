package com.hartwig.hmftools.patientreporter;

import java.io.IOException;

import javax.xml.stream.XMLStreamException;

import com.hartwig.hmftools.patientreporter.algo.AnalysedPatientReport;
import com.hartwig.hmftools.patientreporter.panel.PanelFailReport;
import com.hartwig.hmftools.patientreporter.panel.PanelReport;
import com.hartwig.hmftools.patientreporter.qcfail.QCFailReport;

import org.jetbrains.annotations.NotNull;

public interface ReportWriter {

    void writeAnalysedPatientReport(@NotNull AnalysedPatientReport report, @NotNull String outputFilePath) throws IOException;

    void writeQCFailReport(@NotNull QCFailReport report, @NotNull String outputFilePath) throws IOException;

    void writeJsonFailedFile(@NotNull QCFailReport report, @NotNull String outputFilePath) throws IOException;

    void writeJsonAnalysedFile(@NotNull AnalysedPatientReport report, @NotNull String outputFilePath) throws IOException;

    void writeXMLAnalysedFile(@NotNull AnalysedPatientReport report, @NotNull String outputFilePath) throws IOException, XMLStreamException;

    void writePanelAnalysedReport(@NotNull PanelReport report, @NotNull String outputFilePath) throws IOException;

    void writePanelQCFailReport(@NotNull PanelFailReport report, @NotNull String outputFilePath) throws IOException;

    void writeJsonPanelFile(@NotNull PanelReport report, @NotNull String outputFilePath) throws IOException;

    void writeJsonPanelFailedFile(@NotNull PanelFailReport report, @NotNull String outputFilePath) throws IOException;
}