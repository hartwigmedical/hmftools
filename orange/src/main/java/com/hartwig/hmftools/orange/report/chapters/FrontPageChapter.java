package com.hartwig.hmftools.orange.report.chapters;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;
import static net.sf.dynamicreports.report.builder.DynamicReports.stl;

import com.hartwig.hmftools.orange.report.OrangeColors;
import com.hartwig.hmftools.orange.report.OrangeStyles;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.pdfdata.FrontPageData;
import com.hartwig.hmftools.orange.report.tables.FrontPageTables;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.HorizontalListBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.constant.HorizontalImageAlignment;
import net.sf.dynamicreports.report.constant.PageOrientation;
import net.sf.dynamicreports.report.constant.PageType;
import net.sf.jasperreports.engine.JREmptyDataSource;

public class FrontPageChapter implements ReportChapter
{
    private final FrontPageData mData;

    public FrontPageChapter(final FrontPageData data, final Object unused)
    {
        mData = data;
    }

    @Override
    public String name()
    {
        return "Front Page";
    }

    @Override
    public boolean isLandscape()
    {
        return false;
    }

    @Override
    public JasperReportBuilder buildReport()
    {
        JasperReportBuilder report = report().setPageFormat(PageType.A4, PageOrientation.PORTRAIT);

        JasperReportBuilder sampleTable = FrontPageTables.buildSampleSummary(mData.sampleSummary(), mData.qcWarning());
        JasperReportBuilder pipelineTable = FrontPageTables.buildTechnicalSummary(mData.technicalSummary());

        int rows = Math.max(mData.driverSummary().size(), mData.genomeWideFeatures().size());
        VerticalListBuilder leftBox = FrontPageTables.buildDriverSummaryBox(mData.driverSummary(), rows);
        VerticalListBuilder rightBox = FrontPageTables.buildGenomeWideFeaturesBox(mData.genomeWideFeatures(), rows);
        HorizontalListBuilder tablePair = cmp.horizontalList(leftBox, rightBox);

        VerticalListBuilder tablesContainer = cmp.verticalList(tablePair);
        tablesContainer.setStyle(OrangeStyles.BORDERED_BOX_STYLE);

        var circosImage = cmp.image(mData.circosPlotPath())
                .setHeight(ReportResources.FRONT_CIRCOS_IMAGE_HEIGHT)
                .setHorizontalImageAlignment(HorizontalImageAlignment.CENTER);
        circosImage.setStyle(OrangeStyles.BORDERED_BOX_STYLE);

        VerticalListBuilder componentStack = cmp.verticalList(
                cmp.subreport(sampleTable),
                cmp.line().setPen(stl.penThin().setLineColor(OrangeColors.PALETTE_LIGHT_GREY)),
                cmp.verticalGap(10),
                cmp.subreport(pipelineTable),
                cmp.line().setPen(stl.penThin().setLineColor(OrangeColors.PALETTE_LIGHT_GREY)),
                cmp.verticalGap(10),
                tablesContainer,
                circosImage
        );

        return report.title(componentStack);
    }
}
