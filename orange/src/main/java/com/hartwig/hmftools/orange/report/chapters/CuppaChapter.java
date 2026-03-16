package com.hartwig.hmftools.orange.report.chapters;

import static com.hartwig.hmftools.orange.report.ReportResources.FULL_PAGE_IMAGE_HEIGHT;
import static com.hartwig.hmftools.orange.report.ReportResources.FULL_PAGE_IMAGE_WIDTH;

import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.orange.report.PlotPathResolver;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.datamodel.purple.PurpleQCInterpretation;
import com.hartwig.hmftools.orange.report.util.Images;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.property.HorizontalAlignment;

import org.jetbrains.annotations.NotNull;

public class CuppaChapter implements ReportChapter
{
    private final OrangeRecord mReport;
    private final PlotPathResolver mPlotPathResolver;
    private final ReportResources mReportResources;

    public CuppaChapter(final OrangeRecord report, final PlotPathResolver plotPathResolver, final ReportResources reportResources)
    {
        mReport = report;
        mPlotPathResolver = plotPathResolver;
        mReportResources = reportResources;
    }

    @Override
    public String name()
    {
        return "Tissue of Origin";
    }

    @Override
    public PageSize pageSize()
    {
        return PageSize.A4.rotate();
    }

    @Override
    public void render(final Document document)
    {
        document.add(new Paragraph(name()).addStyle(mReportResources.chapterTitleStyle()));

        boolean isFail = PurpleQCInterpretation.isFail(mReport.purple().fit().qc());
        if(!isFail && mReport.plots().cuppaSummaryPlot() != null)
        {
            addCuppaSummaryPlot(document);
        }
        else
        {
            document.add(new Paragraph(ReportResources.NOT_AVAILABLE).addStyle(mReportResources.tableContentStyle()));
        }
    }

    private void addCuppaSummaryPlot(@NotNull Document document)
    {
        Image cuppaSummaryPlot = Images.build(mPlotPathResolver.resolve(mReport.plots().cuppaSummaryPlot()));
        cuppaSummaryPlot.setMaxWidth(FULL_PAGE_IMAGE_WIDTH);
        cuppaSummaryPlot.setMaxHeight(FULL_PAGE_IMAGE_HEIGHT);
        cuppaSummaryPlot.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(cuppaSummaryPlot);
    }
}
