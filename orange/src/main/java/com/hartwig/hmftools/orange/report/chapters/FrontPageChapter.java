package com.hartwig.hmftools.orange.report.chapters;

import static com.hartwig.hmftools.orange.report.ReportResources.FRONT_CIRCOS_IMAGE_HEIGHT;
import static com.hartwig.hmftools.orange.report.ReportResources.PAGE_MARGIN_LEFT;
import static com.hartwig.hmftools.orange.report.ReportResources.PAGE_MARGIN_RIGHT;

import java.io.IOException;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.OrangeConfig;
import com.hartwig.hmftools.orange.report.DocumentContext;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.tables.FrontPageTables;

import org.apache.pdfbox.pdmodel.PDPageContentStream;
import org.apache.pdfbox.pdmodel.common.PDRectangle;
import org.jetbrains.annotations.NotNull;

import be.quodlibet.boxable.BaseTable;

public class FrontPageChapter implements ReportChapter
{
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

    @NotNull
    @Override
    public String name()
    {
        return "Front Page";
    }

    @NotNull
    @Override
    public PDRectangle pageSize()
    {
        return PDRectangle.A4;
    }

    @Override
    public void render(@NotNull final DocumentContext document) throws IOException
    {
        final float contentWidth = pageSize().getWidth() - (PAGE_MARGIN_LEFT + PAGE_MARGIN_RIGHT);//contentWidth();
        document.addTable(FrontPageTables.buildSampleSummary(document, mReport, mConfig, contentWidth, mReportResources));
        document.addTable(FrontPageTables.buildTechnicalSummary(document, mReport, mConfig, contentWidth, mReportResources));
        document.addSpacing(10);

        // Driver Summary and Genome Wide Biomarkers side-by-side
        float nestedTableWidth = (contentWidth / 2) - 3 * 5;

        float savedY = document.cursorY();
        float driverTableLeft = PAGE_MARGIN_LEFT + 5 + 5;
        BaseTable driverTable = FrontPageTables.buildDriverSummary(document, mReport, nestedTableWidth, driverTableLeft, mReportResources);
        float leftFinalY = document.addTableNoAdvance(driverTable);

        document.setCursorY(savedY); // reset cursor to same line
        float genomeTableLeft = document.marginLeft() + nestedTableWidth + 4 * 5;
        BaseTable genomeTable =
                FrontPageTables.buildGenomeWideFeaturesAtX(document, mReport, nestedTableWidth, genomeTableLeft, mReportResources);
        float rightFinalY = document.addTableNoAdvance(genomeTable);

        // Draw box borders around both sections
        float bottomY = Math.min(leftFinalY, rightFinalY);
        float borderBoxesWidths = nestedTableWidth + 10;
        float borderBoxesHeight = savedY - bottomY;
        drawBoxBorder(document, driverTableLeft - 5, bottomY, borderBoxesWidths, borderBoxesHeight);
        drawBoxBorder(document, genomeTableLeft - 5, bottomY, borderBoxesWidths, borderBoxesHeight);
        drawBoxBorder(document, driverTableLeft - 10, savedY + 5, contentWidth, -1 * FRONT_CIRCOS_IMAGE_HEIGHT - 15 - borderBoxesHeight);

        // Advance cursor to the lower of the two tables
        document.setCursorY(bottomY - 5);

        // Add circos plot
        document.addImage(
                mPlotPathResolver.resolve(mReport.plots().purpleFinalCircosPlot()),
                contentWidth, FRONT_CIRCOS_IMAGE_HEIGHT);
    }

    private static void drawBoxBorder(final DocumentContext document, float x, float y, float width, float height) throws IOException
    {
        try(PDPageContentStream cs = new PDPageContentStream(
                document.document(), document.currentPage(), PDPageContentStream.AppendMode.APPEND, true, true))
        {
            cs.setStrokingColor(ReportResources.PALETTE_BLACK);
            cs.setLineWidth(0.5f);
            cs.addRect(x, y, width, height);
            cs.stroke();
        }
    }
}
