package com.hartwig.hmftools.orange.report.chapters;

import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.orange.report.DocumentContext;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.algo.QcStatusInterpretation;
import com.hartwig.hmftools.orange.report.tables.ExpressionTable;
import com.hartwig.hmftools.orange.report.tables.NovelSpliceJunctionTable;
import com.hartwig.hmftools.orange.report.tables.RnaFusionTable;
import com.hartwig.hmftools.orange.report.tables.RnaStatisticsTable;

import org.apache.pdfbox.pdmodel.common.PDRectangle;
import org.jetbrains.annotations.NotNull;

public class RnaFindingsChapter implements ReportChapter
{
    private final IsofoxRecord mIsofoxRecord;
    private final ReportResources mReportResources;

    public RnaFindingsChapter(final IsofoxRecord isofox, final ReportResources reportResources)
    {
        mIsofoxRecord = isofox;
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
    public PDRectangle pageSize()
    {
        return PDRectangle.A4;
    }

    @Override
    public void render(@NotNull final DocumentContext document) throws IOException
    {
        document.addParagraph(name(), mReportResources.chapterTitleStyle());

        if(QcStatusInterpretation.hasRnaFail(mIsofoxRecord))
        {
            document.addQcFailNotice(mReportResources);
            return;
        }

        addStatistics(document);
        addExpressionTables(document);
        addRnaFusionTables(document);
        addNovelSpliceJunctionTables(document);
    }

    private void addStatistics(final DocumentContext document) throws IOException
    {
        String title = "QC";
        document.addTable(RnaStatisticsTable.build(document, title, contentWidth(), mIsofoxRecord.summary(), mReportResources));
    }

    private void addExpressionTables(final DocumentContext document) throws IOException
    {
        String highExpressionTitle = "High Gene Expression";

        List<GeneExpression> reportableHighExpression = mIsofoxRecord.highExpressionGenes();
        String titleHighExpression = highExpressionTitle + " (" + reportableHighExpression.size() + ")";

        document.addTable(ExpressionTable.build(
                document, titleHighExpression, contentWidth(), reportableHighExpression, false, mReportResources));
    }

    private void addRnaFusionTables(final DocumentContext document) throws IOException
    {
        String fusionsTitle = "Fusions Detected In RNA, Not In DNA";

        List<RnaFusion> reportableNovelKnownFusions = mIsofoxRecord.fusions();
        String titleKnownFusions = fusionsTitle + " (" + reportableNovelKnownFusions.size() + ")";
        document.addTable(RnaFusionTable.build(document, titleKnownFusions, contentWidth(), reportableNovelKnownFusions, mReportResources));
    }

    private void addNovelSpliceJunctionTables(final DocumentContext document) throws IOException
    {
        String novelSplicJunctionsTitle = "Novel Splice Junctions";

        List<NovelSpliceJunction> reportableSkippedExons = mIsofoxRecord.novelSpliceJunctions();
        String titleSkippedExonJunctions = novelSplicJunctionsTitle + " (" + reportableSkippedExons.size() + ")";
        document.addTable(NovelSpliceJunctionTable.build(document, titleSkippedExonJunctions, contentWidth(), reportableSkippedExons, mReportResources));
    }
}
