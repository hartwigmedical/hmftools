package com.hartwig.hmftools.orange.report.chapters;

import static com.hartwig.hmftools.orange.report.ReportResources.FULL_PAGE_IMAGE_HEIGHT;
import static com.hartwig.hmftools.orange.report.ReportResources.FULL_PAGE_IMAGE_WIDTH;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.report.DocumentContext;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;

import org.apache.pdfbox.pdmodel.common.PDRectangle;
import org.jetbrains.annotations.NotNull;

public class QualityControlChapter implements ReportChapter
{
    private static final PDRectangle A4_LANDSCAPE = new PDRectangle(PDRectangle.A4.getHeight(), PDRectangle.A4.getWidth());

    private final OrangeRecord mReport;
    private final PlotPathResolver mPlotPathResolver;
    private final ReportResources mReportResources;

    public QualityControlChapter(final OrangeRecord report, final PlotPathResolver plotPathResolver, final ReportResources reportResources)
    {
        mReport = report;
        mPlotPathResolver = plotPathResolver;
        mReportResources = reportResources;
    }

    @NotNull
    @Override
    public String name()
    {
        return "Quality Control";
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

        List<String> qSeePlots = mReport.plots().qSeePlots();

        if(qSeePlots.isEmpty())
        {
            document.addParagraph("Plot not available", mReportResources.tableContentStyle());
            return;
        }

        for(int i = 0; i < qSeePlots.size(); ++i)
        {
            if(i > 0)
            {
                document.newPage();
                document.addParagraph(name(), mReportResources.chapterTitleStyle());
            }

            document.addImage(mPlotPathResolver.resolve(qSeePlots.get(i)), FULL_PAGE_IMAGE_WIDTH, FULL_PAGE_IMAGE_HEIGHT);
        }
    }
}
