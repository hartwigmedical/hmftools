package com.hartwig.hmftools.orange.report.chapters;

import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;

import org.jetbrains.annotations.NotNull;

public class QualityControlChapter implements ReportChapter {

    @NotNull
    private final OrangeReport report;

    public QualityControlChapter(@NotNull final OrangeReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String name() {
        return "Quality Control";
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph("Quality Control"));
        document.add(new Paragraph("TODO: Add Purple Plots").addStyle(ReportResources.tableContentStyle()));
        document.add(new Paragraph("TODO: Add BQR plots").addStyle(ReportResources.tableContentStyle()));
        document.add(new Paragraph("TODO: Add Metrics / Flagstats").addStyle(ReportResources.tableContentStyle()));
        document.add(new Paragraph("TODO: Add Contamination").addStyle(ReportResources.tableContentStyle()));
    }
}
