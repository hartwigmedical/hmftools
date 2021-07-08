package com.hartwig.hmftools.orange.report.chapters;

import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;

import org.jetbrains.annotations.NotNull;

public class DatabaseCompareChapter implements ReportChapter {

    @NotNull
    private final OrangeReport report;

    public DatabaseCompareChapter(@NotNull final OrangeReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String name() {
        return "Database Compare";
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph("Database compare"));
        document.add(new Paragraph("TODO: Add Charles Cuppa table").addStyle(ReportResources.tableContentStyle()));
    }
}
