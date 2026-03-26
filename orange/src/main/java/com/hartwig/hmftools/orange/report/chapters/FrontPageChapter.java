package com.hartwig.hmftools.orange.report.chapters;

import static java.lang.Math.floor;
import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.orange.report.ReportResources.FRONT_CIRCOS_IMAGE_HEIGHT;
import static com.hartwig.hmftools.orange.report.ReportResources.PAGE_MARGIN_BOTTOM;
import static com.hartwig.hmftools.orange.report.ReportResources.PAGE_MARGIN_TOP;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.tables.FrontPageTables;
import com.hartwig.hmftools.orange.report.util.Images;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.kernel.geom.Rectangle;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.borders.Border;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.layout.LayoutArea;
import com.itextpdf.layout.layout.LayoutContext;
import com.itextpdf.layout.layout.LayoutResult;
import com.itextpdf.layout.property.HorizontalAlignment;
import com.itextpdf.layout.property.UnitValue;
import com.itextpdf.layout.renderer.IRenderer;

public class FrontPageChapter implements ReportChapter
{
    private static final String NONE = "None";

    private final OrangeConfig mConfig;
    private final OrangeRecord mReport;
    private final PlotPathResolver mPlotPathResolver;
    private final ReportResources mReportResources;

    public FrontPageChapter(
            final OrangeConfig config, final OrangeRecord report, final PlotPathResolver plotPathResolver,
            final ReportResources reportResources)
    {
        mConfig = config;
        mReport = report;
        mPlotPathResolver = plotPathResolver;
        mReportResources = reportResources;
    }

    @Override
    public String name()
    {
        return "Front Page";
    }

    @Override
    public PageSize pageSize()
    {
        return PageSize.A4;
    }

    @Override
    public void render(final Document document)
    {
        Table sampleSummaryTable = FrontPageTables.buildSampleSummary(mReport, mConfig, contentWidth(), mReportResources);
        document.add(sampleSummaryTable);

        Table technicalSummaryTable = FrontPageTables.buildTechnicalSummary(mReport, mConfig, contentWidth(), mReportResources);
        technicalSummaryTable.setMarginBottom(10);
        document.add(technicalSummaryTable);

        Table topTable = new Table(UnitValue.createPercentArray(new float[] { 1, 1 })).setWidth(contentWidth() - 5);

        Table driverSummaryTable = FrontPageTables.buildDriverSummary(mReport, contentWidth(), mReportResources);
        Table genomeWideTable = FrontPageTables.buildGenomeWideFeatures(mReport, contentWidth(), mReportResources);

        topTable.addCell(driverSummaryTable);
        topTable.addCell(genomeWideTable);

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1 })).setWidth(contentWidth()).setPadding(0);
        table.addCell(topTable);

        float pageHeight = contentHeight();
        IRenderer renderer = table.createRendererSubTree().setParent(document.getRenderer());
        LayoutResult result = renderer.layout(new LayoutContext(new LayoutArea(0, new Rectangle(contentWidth(), pageHeight))));
        float currentHeight = result.getOccupiedArea().getBBox().getHeight();

        int remainingHeight = (int)floor(pageHeight - currentHeight - PAGE_MARGIN_TOP - PAGE_MARGIN_BOTTOM) - 50;
        int maxCircosHeight = min(remainingHeight, FRONT_CIRCOS_IMAGE_HEIGHT);

        Image circosImage = Images.build(mPlotPathResolver.resolve(mReport.plots().purpleFinalCircosPlot()));
        circosImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        circosImage.setMaxHeight(maxCircosHeight);
        table.addCell(circosImage);

        document.add(table);
    }
}
