package com.hartwig.hmftools.orange.report.chapters;

import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;

import org.jetbrains.annotations.NotNull;

public class CohortComparisonChapter implements ReportChapter {

    @NotNull
    private final OrangeReport report;

    public CohortComparisonChapter(@NotNull final OrangeReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String name() {
        return "Cohort Comparison";
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph("Cohort Comparison").addStyle(ReportResources.chapterTitleStyle()));
        document.add(new Paragraph("TODO: Add cohort predictions for genomic position, signature, driver, expression and alternate splicing")
                .addStyle(ReportResources.tableTitleStyle()));
        document.add(new Paragraph("TODO: Add detailed cohort incidence per driver").addStyle(ReportResources.tableTitleStyle()));
        document.add(new Paragraph("TODO: Add detailed cohort incidence per signature").addStyle(ReportResources.tableTitleStyle()));
    }
}
