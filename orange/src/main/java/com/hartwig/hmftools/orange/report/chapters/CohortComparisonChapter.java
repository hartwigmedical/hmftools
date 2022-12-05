package com.hartwig.hmftools.orange.report.chapters;

import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Images;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.property.HorizontalAlignment;

import org.jetbrains.annotations.NotNull;

public class CohortComparisonChapter implements ReportChapter {

    @NotNull
    private final OrangeReport report;
    @NotNull
    private final PlotPathResolver plotPathResolver;

    public CohortComparisonChapter(@NotNull final OrangeReport report, @NotNull final PlotPathResolver plotPathResolver) {
        this.report = report;
        this.plotPathResolver = plotPathResolver;
    }

    @NotNull
    @Override
    public String name() {
        return "Cohort Comparison";
    }

    @NotNull
    @Override
    public PageSize pageSize() {
        return PageSize.A4.rotate();
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph(name()).addStyle(ReportResources.chapterTitleStyle()));

        addCuppaSummaryPlot(document);
        addCuppaFeaturePlot(document);
    }

    private void addCuppaSummaryPlot(@NotNull Document document) {
        Image cuppaSummaryImage = Images.build(plotPathResolver.resolve(report.plots().cuppaSummaryPlot()));
        cuppaSummaryImage.setMaxWidth(740);
        cuppaSummaryImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(cuppaSummaryImage);
    }

    private void addCuppaFeaturePlot(@NotNull Document document) {
        String featurePlotPaths = report.plots().cuppaFeaturePlot();
        if (featurePlotPaths != null) {
            Image cuppaFeatureImage = Images.build(plotPathResolver.resolve(featurePlotPaths));
            cuppaFeatureImage.setMaxWidth(740);
            cuppaFeatureImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
            document.add(cuppaFeatureImage);
        }
    }
}
