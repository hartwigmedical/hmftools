package com.hartwig.hmftools.orange.report.chapters;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;

import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.pdfdata.QualityControlData;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.constant.HorizontalImageAlignment;
import net.sf.dynamicreports.report.constant.PageOrientation;
import net.sf.dynamicreports.report.constant.PageType;

public class QualityControlChapter implements ReportChapter
{
    private final QualityControlData mData;

    public QualityControlChapter(final QualityControlData data, final Object unused)
    {
        mData = data;
    }

    @Override
    public String name()
    {
        return "Quality Control";
    }

    @Override
    public boolean isLandscape()
    {
        return true;
    }

    @Override
    public JasperReportBuilder buildReport()
    {
        JasperReportBuilder report = report().setPageFormat(PageType.A4, PageOrientation.LANDSCAPE);
        VerticalListBuilder content = cmp.verticalList();
        content.add(cmp.text(name()).setStyle(OrangeFonts.CHAPTER_TITLE_STYLE));

        if(mData.hasPurpleFail)
        {
            content.add(cmp.text(ReportResources.NOT_AVAILABLE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE));
            return report.summary(content);
        }

        if(mData.qSeePlotPath == null)
        {
            content.add(cmp.text("Plot not available").setStyle(OrangeFonts.TABLE_CONTENT_STYLE));
            return report.summary(content);
        }

        content.add(
                cmp.image(mData.qSeePlotPath)
                        .setFixedHeight(ReportResources.FULL_PAGE_IMAGE_HEIGHT)
                        .setHorizontalImageAlignment(HorizontalImageAlignment.CENTER)
        );

        return report.summary(content);
    }
}
