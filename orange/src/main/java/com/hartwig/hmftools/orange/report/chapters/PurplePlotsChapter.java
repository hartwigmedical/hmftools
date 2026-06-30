package com.hartwig.hmftools.orange.report.chapters;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;

import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.pdfdata.PurplePlotsData;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.HorizontalListBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.constant.HorizontalImageAlignment;
import net.sf.dynamicreports.report.constant.PageOrientation;
import net.sf.dynamicreports.report.constant.PageType;

public class PurplePlotsChapter implements ReportChapter
{
    private static final int PLOT_IMAGE_HEIGHT = 225;

    private final PurplePlotsData mData;

    public PurplePlotsChapter(final PurplePlotsData data, final Object unused)
    {
        mData = data;
    }

    @Override
    public String name()
    {
        return "Purity and Ploidy";
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

        // row 1: input circos | CN PDF | somatic clonality
        HorizontalListBuilder row1 = cmp.horizontalList(
                cmp.image(mData.purpleInputCircosPlotPath)
                        .setFixedHeight(PLOT_IMAGE_HEIGHT)
                        .setHorizontalImageAlignment(HorizontalImageAlignment.CENTER),
                cmp.image(mData.purpleCopyNumberPlotPath)
                        .setFixedHeight(PLOT_IMAGE_HEIGHT)
                        .setHorizontalImageAlignment(HorizontalImageAlignment.CENTER),
                cmp.image(mData.purpleClonalityPlotPath)
                        .setFixedHeight(PLOT_IMAGE_HEIGHT)
                        .setHorizontalImageAlignment(HorizontalImageAlignment.CENTER)
        );

        // row 2: ploidy/purity range | minor allele PDF | somatic rainfall
        HorizontalListBuilder row2 = cmp.horizontalList(
                cmp.image(mData.purplePurityRangePlotPath)
                        .setFixedHeight(PLOT_IMAGE_HEIGHT)
                        .setHorizontalImageAlignment(HorizontalImageAlignment.CENTER),
                cmp.image(mData.purpleMinorAlleleMapPlotPath)
                        .setFixedHeight(PLOT_IMAGE_HEIGHT)
                        .setHorizontalImageAlignment(HorizontalImageAlignment.CENTER),
                cmp.image(mData.purpleRainfallPlotPath)
                        .setFixedHeight(PLOT_IMAGE_HEIGHT)
                        .setHorizontalImageAlignment(HorizontalImageAlignment.CENTER)
        );

        content.add(row1);
        content.add(row2);

        return report.summary(content);
    }
}
