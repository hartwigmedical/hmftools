package com.hartwig.hmftools.orange.report.chapters;

import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.ImageUtil;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.property.HorizontalAlignment;

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

        addCuppaReportPlot(document);

    }

    private void addCuppaReportPlot(@NotNull Document document) {
        Image cuppaReportImage = ImageUtil.build(report.plots().cuppaReportPlot());
        cuppaReportImage.setMaxWidth(ReportResources.CONTENT_WIDTH);
        cuppaReportImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(cuppaReportImage);
    }
}
