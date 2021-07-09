package com.hartwig.hmftools.orange.report.chapters;

import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;

import org.jetbrains.annotations.NotNull;

public class GermlineFindingsChapter implements ReportChapter {

    @NotNull
    private final OrangeReport report;
    private final boolean reportGermline;

    public GermlineFindingsChapter(@NotNull final OrangeReport report, final boolean reportGermline) {
        this.report = report;
        this.reportGermline = reportGermline;
    }

    @NotNull
    @Override
    public String name() {
        return "Germline Findings";
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph("Germline Findings"));
        if (reportGermline) {
            document.add(new Paragraph("TODO: Add Germline variants").addStyle(ReportResources.tableContentStyle()));
            document.add(new Paragraph("TODO: Add MVLH Analysis").addStyle(ReportResources.tableContentStyle()));
            document.add(new Paragraph("TODO: Add Germline CN aberrations").addStyle(ReportResources.tableContentStyle()));
            document.add(new Paragraph("TODO: Add PEACH").addStyle(ReportResources.tableContentStyle()));
        } else {
            document.add(new Paragraph("N/A").addStyle(ReportResources.tableContentStyle()));
        }
    }
}
