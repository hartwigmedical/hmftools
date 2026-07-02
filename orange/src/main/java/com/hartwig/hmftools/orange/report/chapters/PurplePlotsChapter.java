package com.hartwig.hmftools.orange.report.chapters;

import java.io.IOException;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.algo.QcStatusInterpretation;
import com.hartwig.hmftools.orange.report.DocumentContext;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.apache.pdfbox.pdmodel.common.PDRectangle;
import org.jetbrains.annotations.NotNull;

public class PurplePlotsChapter implements ReportChapter
{
    private static final PDRectangle A4_LANDSCAPE = new PDRectangle(PDRectangle.A4.getHeight(), PDRectangle.A4.getWidth());
    private static final int PLOT_IMAGE_HEIGHT = 225;

    private final OrangeRecord mReport;
    private final PlotPathResolver mPlotPathResolver;
    private final ReportResources mReportResources;

    public PurplePlotsChapter(final OrangeRecord report, final PlotPathResolver plotPathResolver, final ReportResources reportResources)
    {
        mReport = report;
        mPlotPathResolver = plotPathResolver;
        mReportResources = reportResources;
    }

    @NotNull
    @Override
    public String name()
    {
        return "Purity and Ploidy";
    }

    @NotNull
    @Override
    public PDRectangle pageSize()
    {
        return A4_LANDSCAPE;
    }

    @Override
    public void render(@NotNull final DocumentContext document) throws IOException
    {
        document.addParagraph(name(), mReportResources.chapterTitleStyle());

        if(QcStatusInterpretation.hasPurpleFail(mReport.purple().fit().qc()))
        {
            document.addQcFailNotice(mReportResources);
            return;
        }

        addPurplePlots(document);
    }

    private void addPurplePlots(final DocumentContext document) throws IOException
    {
        // layout: 3 columns x 2 rows on a single landscape page
        // row 1: input circos              CN plot                 somatic clonality
        // row 2: ploidy/purity range       minor allele plot       somatic rainfall

        float colWidth = contentWidth() / 3f;
        float marginLeft = document.marginLeft();

        // Row 1
        float rowY = document.cursorY();
        addPlotAt(document, mReport.plots().purpleInputCircosPlot(), marginLeft, rowY, colWidth);
        addPlotAt(document, mReport.plots().purpleCopyNumberPlot(), marginLeft + colWidth, rowY, colWidth);
        addPlotAt(document, mReport.plots().purpleClonalityPlot(), marginLeft + 2 * colWidth, rowY, colWidth);

        // Row 2
        rowY -= PLOT_IMAGE_HEIGHT + 5;
        addPlotAt(document, mReport.plots().purplePurityRangePlot(), marginLeft, rowY, colWidth);
        addPlotAt(document, mReport.plots().purpleMinorAlleleMapPlot(), marginLeft + colWidth, rowY, colWidth);
        addPlotAt(document, mReport.plots().purpleRainfallPlot(), marginLeft + 2 * colWidth, rowY, colWidth);

        document.setCursorY(rowY - PLOT_IMAGE_HEIGHT - 5);
    }

    private void addPlotAt(final DocumentContext document, final String plotFilename, float x, float y, float width) throws IOException
    {
        document.addImageAt(mPlotPathResolver.resolve(plotFilename), x, y, width, PLOT_IMAGE_HEIGHT);
    }
}
