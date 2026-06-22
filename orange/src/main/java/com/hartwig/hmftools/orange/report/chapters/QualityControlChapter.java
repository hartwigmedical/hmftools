package com.hartwig.hmftools.orange.report.chapters;

import static com.hartwig.hmftools.orange.report.ReportResources.FULL_PAGE_IMAGE_HEIGHT;
import static com.hartwig.hmftools.orange.report.ReportResources.FULL_PAGE_IMAGE_WIDTH;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.pdfdata.QualityControlData;
import com.hartwig.hmftools.orange.report.util.Images;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Image;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.property.HorizontalAlignment;

import org.jetbrains.annotations.NotNull;

public class QualityControlChapter implements ReportChapter
{
    private final QualityControlData mData;
    private final ReportResources mReportResources;

    public QualityControlChapter(final QualityControlData data, final ReportResources reportResources)
    {
        mData = data;
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

        if(mData.qSeePlotPath == null)
        {
            document.add(new Paragraph("Plot not available").addStyle(mReportResources.tableContentStyle()));
            return;
        }

        addQSeePdf(document);
    }

    private void addQSeePdf(final Document document)
    {
        Image qSeePdf = Images.build(mData.qSeePlotPath);
        qSeePdf.setMaxWidth(FULL_PAGE_IMAGE_WIDTH);
        qSeePdf.setMaxHeight(FULL_PAGE_IMAGE_HEIGHT);
        qSeePdf.setHorizontalAlignment(HorizontalAlignment.CENTER);
        document.add(qSeePdf);
    }
}
