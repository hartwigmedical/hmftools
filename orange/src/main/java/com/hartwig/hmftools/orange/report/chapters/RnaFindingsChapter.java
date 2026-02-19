package com.hartwig.hmftools.orange.report.chapters;

import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentage;

import java.util.List;
import java.util.Set;
import java.util.StringJoiner;

import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.datamodel.isofox.RnaQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.datamodel.purple.PurpleQCInterpretation;
import com.hartwig.hmftools.orange.report.tables.ExpressionTable;
import com.hartwig.hmftools.orange.report.tables.NovelSpliceJunctionTable;
import com.hartwig.hmftools.orange.report.tables.RnaFusionTable;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

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
        Cells cells = new Cells(mReportResources);
        Table table = Tables.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1 },
                new Cell[] { cells.createHeader("QC"), cells.createHeader("Total Fragments"), cells.createHeader("Non-Duplicate Fragments"),
                        cells.createHeader("Duplicate rate") });

        if(PurpleQCInterpretation.isContaminated(mPurpleRecord.fit().qc()))
        {
            table.addCell(cells.createSpanningEntry(table, ReportResources.NOT_AVAILABLE));
        }
        else
        {
            table.addCell(cells.createContent(qcValue(mIsofoxRecord.summary().qcStatus())));
            table.addCell(cells.createContent(String.valueOf(mIsofoxRecord.summary().totalFragments())));

            long nonDuplicates = mIsofoxRecord.summary().totalFragments() - mIsofoxRecord.summary().duplicateFragments();
            table.addCell(cells.createContent(String.valueOf(nonDuplicates)));

            double duplicateRate = mIsofoxRecord.summary().duplicateFragments() / (double) mIsofoxRecord.summary().totalFragments();
            table.addCell(cells.createContent(formatPercentage(duplicateRate)));

            addQCWarningInCaseOfFail(table, cells);
        }

        document.add(new Tables(mReportResources).createWrapping(table));
    }

    @NotNull
    private static String qcValue(Set<RnaQCStatus> qcStatus)
    {
        StringJoiner joiner = new StringJoiner(", ");
        for(RnaQCStatus status : qcStatus)
        {
            joiner.add(status.name());
        }
        return joiner.toString();
    }

    private void addQCWarningInCaseOfFail(Table table, Cells cells)
    {
        boolean isRnaFail = !mIsofoxRecord.summary().qcStatus().contains(RnaQCStatus.PASS);
        boolean isDnaFailNoTumor = PurpleQCInterpretation.isFailNoTumor(mPurpleRecord.fit().qc());

        if(isRnaFail || isDnaFailNoTumor)
        {
            String warning = isRnaFail ?
                    "The RNA QC status of this sample is not a pass. All presented RNA data should be interpreted with caution"
                    : "The DNA QC status of this sample is fail (no tumor). "
                            + "In addition to DNA findings, all RNA findings should be interpreted with caution";

            table.addCell(cells.createSpanningWarning(table, warning));
        }
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
            List<PurpleGeneCopyNumber> somaticGeneCopyNumbers = mPurpleRecord.somaticGeneCopyNumbers();

            List<GeneExpression> reportableHighExpression = mIsofoxRecord.reportableHighExpression();
            String titleHighExpression = highExpressionTitle + " (" + reportableHighExpression.size() + ")";
            document.add(ExpressionTable.build(titleHighExpression, contentWidth(), reportableHighExpression, false, somaticGeneCopyNumbers,
                    mReportResources));

            List<GeneExpression> reportableLowExpression = mIsofoxRecord.reportableLowExpression();
            String titleLowExpression = lowExpressionTitle + " (" + reportableLowExpression.size() + ")";
            document.add(ExpressionTable.build(titleLowExpression, contentWidth(), reportableLowExpression, true, somaticGeneCopyNumbers,
                    mReportResources));
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
            List<RnaFusion> reportableNovelKnownFusions = mIsofoxRecord.reportableNovelKnownFusions();
            String titleKnownFusions = knownFusionsTitle + " (" + reportableNovelKnownFusions.size() + ")";
            document.add(RnaFusionTable.build(titleKnownFusions, contentWidth(), reportableNovelKnownFusions, mReportResources));

            List<RnaFusion> reportableNovelPromiscuous = mIsofoxRecord.reportableNovelPromiscuousFusions();
            String titlePromiscuousFusions = promiscuousFusionsTitle + " (" + reportableNovelPromiscuous.size() + ")";
            document.add(RnaFusionTable.build(titlePromiscuousFusions, contentWidth(), reportableNovelPromiscuous, mReportResources));
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
            List<NovelSpliceJunction> reportableSkippedExons = mIsofoxRecord.reportableSkippedExons();
            String titleSkippedExonJunctions = skippedExonsTitle + " (" + reportableSkippedExons.size() + ")";
            document.add(NovelSpliceJunctionTable.build(titleSkippedExonJunctions, contentWidth(), reportableSkippedExons, mReportResources));

            List<NovelSpliceJunction> reportableNovelExonsIntrons = mIsofoxRecord.reportableNovelExonsIntrons();
            String titleNovelExonIntronJunctions = novelExonsIntronsTitle + " (" + reportableNovelExonsIntrons.size() + ")";
            document.add(NovelSpliceJunctionTable.build(titleNovelExonIntronJunctions, contentWidth(), reportableNovelExonsIntrons,
                    mReportResources));
        }
    }
}
