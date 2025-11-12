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
import com.hartwig.hmftools.datamodel.interpretation.PurpleQCInterpretation;
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
    @NotNull
    private final IsofoxRecord isofox;
    @NotNull
    private final PurpleRecord purple;
    @NotNull
    private final ReportResources reportResources;

    public RnaFindingsChapter(@NotNull final IsofoxRecord isofox, @NotNull final PurpleRecord purple,
            @NotNull final ReportResources reportResources)
    {
        this.isofox = isofox;
        this.purple = purple;
        this.reportResources = reportResources;
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
        document.add(new Paragraph(name()).addStyle(reportResources.chapterTitleStyle()));

        addKeyQC(document);
        addExpressionTables(document);
        addRnaFusionTables(document);
        addNovelSpliceJunctionTables(document);
    }

    private void addKeyQC(@NotNull Document document)
    {
        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1 },
                new Cell[] { cells.createHeader("QC"), cells.createHeader("Total Fragments"), cells.createHeader("Non-Duplicate Fragments"),
                        cells.createHeader("Duplicate rate") });

        if(PurpleQCInterpretation.isContaminated(purple.fit().qc()))
        {
            table.addCell(cells.createSpanningEntry(table, ReportResources.NOT_AVAILABLE));
        }
        else
        {
            table.addCell(cells.createContent(qcValue(isofox.summary().qcStatus())));
            table.addCell(cells.createContent(String.valueOf(isofox.summary().totalFragments())));

            long nonDuplicates = isofox.summary().totalFragments() - isofox.summary().duplicateFragments();
            table.addCell(cells.createContent(String.valueOf(nonDuplicates)));

            double duplicateRate = isofox.summary().duplicateFragments() / (double) isofox.summary().totalFragments();
            table.addCell(cells.createContent(formatPercentage(duplicateRate)));

            addQCWarningInCaseOfFail(table, cells);
        }

        document.add(new Tables(reportResources).createWrapping(table));
    }

    @NotNull
    private static String qcValue(@NotNull Set<RnaQCStatus> qcStatus)
    {
        StringJoiner joiner = new StringJoiner(", ");
        for(RnaQCStatus status : qcStatus)
        {
            joiner.add(status.name());
        }
        return joiner.toString();
    }

    private void addQCWarningInCaseOfFail(@NotNull Table table, @NotNull Cells cells)
    {
        boolean isRnaFail = !isofox.summary().qcStatus().contains(RnaQCStatus.PASS);
        boolean isDnaFailNoTumor = PurpleQCInterpretation.isFailNoTumor(purple.fit().qc());

        if(isRnaFail || isDnaFailNoTumor)
        {
            String warning = isRnaFail ?
                    "The RNA QC status of this sample is not a pass. All presented RNA data should be interpreted with caution"
                    : "The DNA QC status of this sample is fail (no tumor). "
                            + "In addition to DNA findings, all RNA findings should be interpreted with caution";

            table.addCell(cells.createSpanningWarning(table, warning));
        }
    }

    private void addExpressionTables(@NotNull Document document)
    {
        String highExpressionTitle = "Genes with high expression";
        String lowExpressionTitle = "Genes with low expression";

        if(PurpleQCInterpretation.isContaminated(purple.fit().qc()))
        {
            Tables tables = new Tables(reportResources);
            document.add(tables.createNotAvailable(highExpressionTitle, contentWidth()));
            document.add(tables.createNotAvailable(lowExpressionTitle, contentWidth()));
        }
        else
        {
            List<PurpleGeneCopyNumber> somaticGeneCopyNumbers = purple.allSomaticGeneCopyNumbers();

            List<GeneExpression> reportableHighExpression = isofox.reportableHighExpression();
            String titleHighExpression = highExpressionTitle + " (" + reportableHighExpression.size() + ")";
            document.add(ExpressionTable.build(titleHighExpression, contentWidth(), reportableHighExpression, false, somaticGeneCopyNumbers,
                    reportResources));

            List<GeneExpression> reportableLowExpression = isofox.reportableLowExpression();
            String titleLowExpression = lowExpressionTitle + " (" + reportableLowExpression.size() + ")";
            document.add(ExpressionTable.build(titleLowExpression, contentWidth(), reportableLowExpression, true, somaticGeneCopyNumbers,
                    reportResources));
        }
    }

    private void addRnaFusionTables(@NotNull Document document)
    {
        String knownFusionsTitle = "Known fusions detected in RNA and not in DNA";
        String promiscuousFusionsTitle = "Promiscuous fusions detected in RNA and not in DNA";

        if(PurpleQCInterpretation.isContaminated(purple.fit().qc()))
        {
            Tables tables = new Tables(reportResources);
            document.add(tables.createNotAvailable(knownFusionsTitle, contentWidth()));
            document.add(tables.createNotAvailable(promiscuousFusionsTitle, contentWidth()));
        }
        else
        {
            List<RnaFusion> reportableNovelKnownFusions = isofox.reportableNovelKnownFusions();
            String titleKnownFusions = knownFusionsTitle + " (" + reportableNovelKnownFusions.size() + ")";
            document.add(RnaFusionTable.build(titleKnownFusions, contentWidth(), reportableNovelKnownFusions, reportResources));

            List<RnaFusion> reportableNovelPromiscuous = isofox.reportableNovelPromiscuousFusions();
            String titlePromiscuousFusions = promiscuousFusionsTitle + " (" + reportableNovelPromiscuous.size() + ")";
            document.add(RnaFusionTable.build(titlePromiscuousFusions, contentWidth(), reportableNovelPromiscuous, reportResources));
        }
    }

    private void addNovelSpliceJunctionTables(@NotNull Document document)
    {
        String skippedExonsTitle = "Potentially interesting novel splice junctions - Skipped exons";
        String novelExonsIntronsTitle = "Potentially interesting novel splice junctions - Novel exon/intron";

        if(PurpleQCInterpretation.isContaminated(purple.fit().qc()))
        {
            Tables tables = new Tables(reportResources);
            document.add(tables.createNotAvailable(skippedExonsTitle, contentWidth()));
            document.add(tables.createNotAvailable(novelExonsIntronsTitle, contentWidth()));
        }
        else
        {
            List<NovelSpliceJunction> reportableSkippedExons = isofox.reportableSkippedExons();
            String titleSkippedExonJunctions = skippedExonsTitle + " (" + reportableSkippedExons.size() + ")";
            document.add(NovelSpliceJunctionTable.build(titleSkippedExonJunctions, contentWidth(), reportableSkippedExons, reportResources));

            List<NovelSpliceJunction> reportableNovelExonsIntrons = isofox.reportableNovelExonsIntrons();
            String titleNovelExonIntronJunctions = novelExonsIntronsTitle + " (" + reportableNovelExonsIntrons.size() + ")";
            document.add(NovelSpliceJunctionTable.build(titleNovelExonIntronJunctions, contentWidth(), reportableNovelExonsIntrons,
                    reportResources));
        }
    }
}
