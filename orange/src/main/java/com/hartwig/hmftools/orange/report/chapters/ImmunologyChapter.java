package com.hartwig.hmftools.orange.report.chapters;

import com.hartwig.hmftools.orange.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;

import org.jetbrains.annotations.NotNull;

public class ImmunologyChapter implements ReportChapter {

    @NotNull
    private final OrangeReport report;

    public ImmunologyChapter(@NotNull final OrangeReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String name() {
        return "Immunology";
    }

    @NotNull
    @Override
    public PageSize pageSize() {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph(name()).addStyle(ReportResources.chapterTitleStyle()));
        document.add(new Paragraph("Coming soon: HLA types for patient + status in tumor").addStyle(ReportResources.tableTitleStyle()));
        document.add(new Paragraph("Coming soon: List of neo-epitopes with predicted binding affinity").addStyle(ReportResources.tableTitleStyle()));
        document.add(new Paragraph("Coming soon: Details about RNA tumor micro-environment").addStyle(ReportResources.tableTitleStyle()));
    }
}
