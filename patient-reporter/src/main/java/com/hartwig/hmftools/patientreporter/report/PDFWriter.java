package com.hartwig.hmftools.patientreporter.report;

import com.hartwig.hmftools.patientreporter.PatientReport;

import org.jetbrains.annotations.NotNull;

public class PDFWriter {

    @NotNull
    private final String outputDirectory;
    @NotNull
    private final String hmfLogo;

    public PDFWriter(@NotNull final String outputDirectory, @NotNull final String hmfLogo) {
        this.outputDirectory = outputDirectory;
        this.hmfLogo = hmfLogo;
    }

    public void write(@NotNull final PatientReport patientReport) {
        // KODU: TODO!
    }
}
