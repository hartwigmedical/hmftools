package com.hartwig.hmftools.patientreporter.cfreport.chapters;

import com.hartwig.hmftools.patientreporter.QCFailReport;
import com.hartwig.hmftools.patientreporter.cfreport.components.TumorLocationAndTypeTable;
import com.itextpdf.layout.Document;
import org.jetbrains.annotations.NotNull;

import java.io.IOException;

public class QCFailChapter implements ReportChapter {

    private final QCFailReport failReport;

    public QCFailChapter(QCFailReport failReport) {
        this.failReport = failReport;
    }

    @Override
    public String getName() {
        return failReport.reason().title();
    }

    public boolean isFullWidth() {
        return false;
    }

    public boolean hasCompleteSidebar() {
        return true;
    }

    @Override
    public void render(@NotNull Document reportDocument) throws IOException {

        reportDocument.add(TumorLocationAndTypeTable.createTumorLocationAndType(
                failReport.sampleReport().primaryTumorLocationString(),
                failReport.sampleReport().cancerSubTypeString(),
                getContentWidth()));

    }

}
