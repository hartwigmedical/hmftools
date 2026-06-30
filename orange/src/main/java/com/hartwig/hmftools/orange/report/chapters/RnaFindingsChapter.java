package com.hartwig.hmftools.orange.report.chapters;

import static net.sf.dynamicreports.report.builder.DynamicReports.cmp;
import static net.sf.dynamicreports.report.builder.DynamicReports.report;

import com.hartwig.hmftools.orange.report.OrangeFonts;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.pdfdata.RnaFindingsData;
import com.hartwig.hmftools.orange.report.tables.ExpressionTable;
import com.hartwig.hmftools.orange.report.tables.NovelSpliceJunctionTable;
import com.hartwig.hmftools.orange.report.tables.RnaFusionTable;
import com.hartwig.hmftools.orange.report.tables.RnaStatisticsTable;

import net.sf.dynamicreports.jasper.builder.JasperReportBuilder;
import net.sf.dynamicreports.report.builder.component.VerticalListBuilder;
import net.sf.dynamicreports.report.constant.PageOrientation;
import net.sf.dynamicreports.report.constant.PageType;

public class RnaFindingsChapter implements ReportChapter
{
    private final RnaFindingsData mData;

    public RnaFindingsChapter(final RnaFindingsData data, final Object unused)
    {
        mData = data;
    }

    @Override
    public String name()
    {
        return "RNA Findings";
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
        VerticalListBuilder content = cmp.verticalList();
        content.add(cmp.text(name()).setStyle(OrangeFonts.CHAPTER_TITLE_STYLE));

        if(mData.hasRnaFail)
        {
            content.add(cmp.text(ReportResources.NOT_AVAILABLE).setStyle(OrangeFonts.TABLE_CONTENT_STYLE));
            return report.summary(content);
        }

        content.add(cmp.subreport(RnaStatisticsTable.build("QC", mData.summary, null)));

        String titleHighExpression = "High Gene Expression (" + mData.highExpressionGenes.size() + ")";
        content.add(cmp.subreport(ExpressionTable.build(titleHighExpression, mData.highExpressionGenes, false, null)));

        String titleFusions = "Fusions Detected In RNA, Not In DNA (" + mData.fusions.size() + ")";
        content.add(cmp.subreport(RnaFusionTable.build(titleFusions, mData.fusions, null)));

        String titleSpliceJunctions = "Novel Splice Junctions (" + mData.novelSpliceJunctions.size() + ")";
        content.add(cmp.subreport(NovelSpliceJunctionTable.build(titleSpliceJunctions, mData.novelSpliceJunctions, null)));

        return report.summary(content);
    }
}
