package com.hartwig.hmftools.orange.report.chapters;

import static com.hartwig.hmftools.orange.report.ReportResources.FULL_PAGE_IMAGE_HEIGHT;
import static com.hartwig.hmftools.orange.report.ReportResources.FULL_PAGE_IMAGE_WIDTH;

import java.io.IOException;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.algo.QcStatusInterpretation;
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

        if(QcStatusInterpretation.hasPurpleFail(mReport.purple().fit().qc()))
        {
            document.addQcFailNotice(mReportResources);
            return;
        }

        if(mReport.plots().qSeePlot() == null)
        {
            document.addParagraph("Plot not available", mReportResources.tableContentStyle());
            return;
        }

        document.addImage(mPlotPathResolver.resolve(mReport.plots().qSeePlot()), FULL_PAGE_IMAGE_WIDTH, FULL_PAGE_IMAGE_HEIGHT);
    }
}
