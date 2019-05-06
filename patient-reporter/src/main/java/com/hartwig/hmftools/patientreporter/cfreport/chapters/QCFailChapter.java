package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.QCFailReport;
import com.itextpdf.layout.Document;
import org.jetbrains.annotations.NotNull;

import java.io.IOException;

public class QCFailChapter implements ReportChapter {

    private final QCFailReport report;

    public QCFailChapter(QCFailReport report) {
        this.report = report;
    }

    @Override
    @NotNull
    public String getName() {
        return "QC failure";
    }

    public boolean isFullWidth() {
        return false;
    }

    public boolean hasCompleteSidebar() {
        return true;
    }

    @Override
    public void render(@NotNull Document reportDocument) {
        //@TODO Use content from this.report
    }

}
