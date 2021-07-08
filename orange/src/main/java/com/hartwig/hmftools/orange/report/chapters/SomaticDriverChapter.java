package com.hartwig.hmftools.orange.report.chapters;

import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;

import org.jetbrains.annotations.NotNull;

public class SomaticDriverChapter implements ReportChapter {

    @NotNull
    private final OrangeReport report;

    public SomaticDriverChapter(@NotNull final OrangeReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String name() {
        return "Somatic Drivers";
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph("Somatic Drivers"));
        document.add(new Paragraph("TODO: Add Somatic variants").addStyle(ReportResources.tableContentStyle()));
        document.add(new Paragraph("TODO: Add Somatic AMPs/DELs").addStyle(ReportResources.tableContentStyle()));
        document.add(new Paragraph("TODO: Add Somatic Disruptions").addStyle(ReportResources.tableContentStyle()));
        document.add(new Paragraph("TODO: Add Somatic Fusions").addStyle(ReportResources.tableContentStyle()));
        document.add(new Paragraph("TODO: Add Somatic Viral Presence").addStyle(ReportResources.tableContentStyle()));
        document.add(new Paragraph("TODO: Add LINX Visualisations").addStyle(ReportResources.tableContentStyle()));
    }
}
