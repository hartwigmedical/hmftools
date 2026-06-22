package com.hartwig.hmftools.orange.report.chapters;

import java.util.List;

import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.pdfdata.RnaFindingsData;
import com.hartwig.hmftools.orange.report.tables.ExpressionTable;
import com.hartwig.hmftools.orange.report.tables.NovelSpliceJunctionTable;
import com.hartwig.hmftools.orange.report.tables.RnaFusionTable;
import com.hartwig.hmftools.orange.report.tables.RnaStatisticsTable;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;

import org.jetbrains.annotations.NotNull;

public class RnaFindingsChapter implements ReportChapter
{
    private final RnaFindingsData mData;
    private final ReportResources mReportResources;

    public RnaFindingsChapter(final RnaFindingsData data, final ReportResources reportResources)
    {
        mData = data;
        mReportResources = reportResources;
    }

    @NotNull
    @Override
    public String name()
    {
        return "RNA Findings";
    }

    @NotNull
    @Override
    public PageSize pageSize()
    {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull final Document document)
    {
        document.add(new Paragraph(name()).addStyle(mReportResources.chapterTitleStyle()));

        if(mData.hasRnaFail)
        {
            mReportResources.addQcFailNotice(document);
            return;
        }

        addStatistics(document);
        addExpressionTables(document);
        addRnaFusionTables(document);
        addNovelSpliceJunctionTables(document);
    }

    private void addStatistics(final Document document)
    {
        document.add(RnaStatisticsTable.build("QC", contentWidth(), mData.summary, mReportResources));
    }

    private void addExpressionTables(final Document document)
    {
        List<GeneExpression> reportableHighExpression = mData.highExpressionGenes;
        String titleHighExpression = "High Gene Expression (" + reportableHighExpression.size() + ")";

        document.add(ExpressionTable.build(
                titleHighExpression, contentWidth(), reportableHighExpression, false, mReportResources));
    }

    private void addRnaFusionTables(final Document document)
    {
        List<RnaFusion> reportableNovelKnownFusions = mData.fusions;
        String titleKnownFusions = "Fusions Detected In RNA, Not In DNA (" + reportableNovelKnownFusions.size() + ")";
        document.add(RnaFusionTable.build(titleKnownFusions, contentWidth(), reportableNovelKnownFusions, mReportResources));
    }

    private void addNovelSpliceJunctionTables(final Document document)
    {
        List<NovelSpliceJunction> reportableSkippedExons = mData.novelSpliceJunctions;
        String titleSkippedExonJunctions = "Novel Splice Junctions (" + reportableSkippedExons.size() + ")";
        document.add(NovelSpliceJunctionTable.build(titleSkippedExonJunctions, contentWidth(), reportableSkippedExons, mReportResources));
    }
}
