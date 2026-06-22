package com.hartwig.hmftools.orange.report.chapters;

import static com.hartwig.hmftools.orange.report.ReportResources.FULL_PAGE_IMAGE_HEIGHT;
import static com.hartwig.hmftools.orange.report.ReportResources.FULL_PAGE_IMAGE_WIDTH;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.pdfdata.CuppaChapterData;
import com.hartwig.hmftools.orange.report.util.Images;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.property.HorizontalAlignment;

import org.jetbrains.annotations.NotNull;

public class CuppaChapter implements ReportChapter
{
    private final CuppaChapterData mData;
    private final ReportResources mReportResources;

    public CuppaChapter(final CuppaChapterData data, final ReportResources reportResources)
    {
        mData = data;
        mReportResources = reportResources;
    }

    @NotNull
    @Override
    public String name()
    {
        return "Tissue of Origin";
    }

    @NotNull
    @Override
    public PageSize pageSize()
    {
        return PageSize.A4.rotate();
    }

    @Override
    public void render(@NotNull final Document document)
    {
        document.add(new Paragraph(name()).addStyle(mReportResources.chapterTitleStyle()));

        if(mData.hasPurpleFail)
        {
            mReportResources.addQcFailNotice(document);
            return;
        }

        if(mData.cuppaSummaryPlotPath == null)
        {
            document.add(new Paragraph(ReportResources.NOT_AVAILABLE).addStyle(mReportResources.tableContentStyle()));
            return;
        }

        addCuppaSummaryPlot(document);
    }

    private void addCuppaSummaryPlot(final Document document)
    {
        Image cuppaSummaryPlot = Images.build(mData.cuppaSummaryPlotPath);
        cuppaSummaryPlot.setMaxWidth(FULL_PAGE_IMAGE_WIDTH);
        cuppaSummaryPlot.setMaxHeight(FULL_PAGE_IMAGE_HEIGHT);
        cuppaSummaryPlot.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(cuppaSummaryPlot);
    }
}
