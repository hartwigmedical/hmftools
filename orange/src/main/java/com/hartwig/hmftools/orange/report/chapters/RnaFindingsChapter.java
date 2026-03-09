package com.hartwig.hmftools.orange.report.chapters;

import static com.hartwig.hmftools.orange.report.tables.ExpressionTable.buildRnaSummary;

import java.util.List;

import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.datamodel.purple.PurpleQCInterpretation;
import com.hartwig.hmftools.orange.report.tables.ExpressionTable;
import com.hartwig.hmftools.orange.report.tables.NovelSpliceJunctionTable;
import com.hartwig.hmftools.orange.report.tables.RnaFusionTable;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Paragraph;

public class RnaFindingsChapter implements ReportChapter
{
    private final IsofoxRecord mIsofoxRecord;
    private final PurpleRecord mPurpleRecord;
    private final ReportResources mReportResources;

    public RnaFindingsChapter(final IsofoxRecord isofox, final PurpleRecord purple, final ReportResources reportResources)
    {
        mIsofoxRecord = isofox;
        mPurpleRecord = purple;
        mReportResources = reportResources;
    }

    @Override
    public String name()
    {
        return "RNA Findings";
    }

    @Override
    public PageSize pageSize()
    {
        return PageSize.A4;
    }

    @Override
    public void render(final Document document)
    {
        document.add(new Paragraph(name()).addStyle(mReportResources.chapterTitleStyle()));

        addKeyQC(document);
        addExpressionTables(document);
        addRnaFusionTables(document);
        addNovelSpliceJunctionTables(document);
    }

    private void addKeyQC(final Document document)
    {
        String title = "QC";

        document.add(buildRnaSummary(title, contentWidth(), mIsofoxRecord.summary(), mReportResources));
    }

    private void addExpressionTables(final Document document)
    {
        String highExpressionTitle = "Genes with high expression";
        String lowExpressionTitle = "Genes with low expression";

        if(PurpleQCInterpretation.isContaminated(mPurpleRecord.fit().qc()))
        {
            Tables tables = new Tables(mReportResources);
            document.add(tables.createNotAvailable(highExpressionTitle, contentWidth()));
            document.add(tables.createNotAvailable(lowExpressionTitle, contentWidth()));
        }
        else
        {
            List<GeneExpression> reportableHighExpression = mIsofoxRecord.highExpressionGenes();
            String titleHighExpression = highExpressionTitle + " (" + reportableHighExpression.size() + ")";

            document.add(ExpressionTable.build(
                    titleHighExpression, contentWidth(), reportableHighExpression, false, mReportResources));

            List<GeneExpression> reportableLowExpression = mIsofoxRecord.lowExpressionGenes();
            String titleLowExpression = lowExpressionTitle + " (" + reportableLowExpression.size() + ")";

            document.add(ExpressionTable.build(
                    titleLowExpression, contentWidth(), reportableLowExpression, true, mReportResources));
        }
    }

    private void addRnaFusionTables(final Document document)
    {
        String knownFusionsTitle = "Known fusions detected in RNA and not in DNA";
        String promiscuousFusionsTitle = "Promiscuous fusions detected in RNA and not in DNA";

        if(PurpleQCInterpretation.isContaminated(mPurpleRecord.fit().qc()))
        {
            Tables tables = new Tables(mReportResources);
            document.add(tables.createNotAvailable(knownFusionsTitle, contentWidth()));
            document.add(tables.createNotAvailable(promiscuousFusionsTitle, contentWidth()));
        }
        else
        {
            List<RnaFusion> reportableNovelKnownFusions = mIsofoxRecord.fusions();
            String titleKnownFusions = knownFusionsTitle + " (" + reportableNovelKnownFusions.size() + ")";
            document.add(RnaFusionTable.build(titleKnownFusions, contentWidth(), reportableNovelKnownFusions, mReportResources));

            /*
            List<RnaFusion> reportableNovelPromiscuous = mIsofoxRecord.reportableNovelPromiscuousFusions();
            String titlePromiscuousFusions = promiscuousFusionsTitle + " (" + reportableNovelPromiscuous.size() + ")";
            document.add(RnaFusionTable.build(titlePromiscuousFusions, contentWidth(), reportableNovelPromiscuous, mReportResources));
            */
        }
    }

    private void addNovelSpliceJunctionTables(final Document document)
    {
        String skippedExonsTitle = "Potentially interesting novel splice junctions - Skipped exons";
        String novelExonsIntronsTitle = "Potentially interesting novel splice junctions - Novel exon/intron";

        if(PurpleQCInterpretation.isContaminated(mPurpleRecord.fit().qc()))
        {
            Tables tables = new Tables(mReportResources);
            document.add(tables.createNotAvailable(skippedExonsTitle, contentWidth()));
            document.add(tables.createNotAvailable(novelExonsIntronsTitle, contentWidth()));
        }
        else
        {
            List<NovelSpliceJunction> reportableSkippedExons = mIsofoxRecord.novelSpliceJunctions();
            String titleSkippedExonJunctions = skippedExonsTitle + " (" + reportableSkippedExons.size() + ")";
            document.add(NovelSpliceJunctionTable.build(titleSkippedExonJunctions, contentWidth(), reportableSkippedExons, mReportResources));

            /*
            List<NovelSpliceJunction> reportableNovelExonsIntrons = mIsofoxRecord.reportableNovelExonsIntrons();
            String titleNovelExonIntronJunctions = novelExonsIntronsTitle + " (" + reportableNovelExonsIntrons.size() + ")";
            document.add(NovelSpliceJunctionTable.build(titleNovelExonIntronJunctions, contentWidth(), reportableNovelExonsIntrons,
                    mReportResources));
            */
        }
    }
}
