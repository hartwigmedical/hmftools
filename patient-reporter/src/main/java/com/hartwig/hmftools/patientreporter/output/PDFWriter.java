package com.hartwig.hmftools.patientreporter.output;

import com.hartwig.hmftools.patientreporter.PatientReport;

import org.jetbrains.annotations.NotNull;

import net.sf.dynamicreports.report.exception.DRException;

public final class PDFWriter {

    private PDFWriter() {
    }

    public static void writeToPDF(@NotNull final String basePath, @NotNull final PatientReport patientReport)
            throws DRException {
        // KODU: TODO!
    }
}
