package com.hartwig.hmftools.orange.report.chapters;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.tables.FrontPageTables;
import com.hartwig.hmftools.orange.report.util.Images;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Table;
import com.itextpdf.layout.property.HorizontalAlignment;
import com.itextpdf.layout.property.UnitValue;

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
        // addSummaryTable(document);
        // addDetailsAndPlots(document);

        Table sampleSummaryTable = FrontPageTables.buildSampleSummary(mReport, mConfig, contentWidth(), mReportResources);

        document.add(new Tables(mReportResources).createWrapping(sampleSummaryTable));

        Table topTable = new Table(UnitValue.createPercentArray(new float[] { 1, 1 })).setWidth(contentWidth() - 5);

        Table driverSummaryTable = FrontPageTables.buildDriverSummary(mReport, contentWidth(), mReportResources);
        Table genomeWideTable = FrontPageTables.buildGenomeWideFeatures(mReport, contentWidth(), mReportResources);

        topTable.addCell(driverSummaryTable);
        topTable.addCell(genomeWideTable);

        Table table = new Table(UnitValue.createPercentArray(new float[] { 1 })).setWidth(contentWidth()).setPadding(0);
        table.addCell(topTable);

        Image circosImage = Images.build(mPlotPathResolver.resolve(mReport.plots().purpleFinalCircosPlot()));
        circosImage.setHorizontalAlignment(HorizontalAlignment.CENTER);
        circosImage.setMaxHeight(400);
        table.addCell(circosImage);

        document.add(table);
    }
}
