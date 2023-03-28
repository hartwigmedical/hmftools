package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.tables.ExpressionTable;
import com.hartwig.hmftools.orange.report.tables.NovelSpliceJunctionTable;
import com.hartwig.hmftools.orange.report.tables.RNAFusionTable;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.kernel.geom.PageSize;
import com.itextpdf.layout.Document;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Paragraph;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public class RNAFindingsChapter implements ReportChapter {

    private static final DecimalFormat PERCENTAGE_FORMAT = ReportResources.decimalFormat("#'%'");

    @NotNull
    private final OrangeRecord report;
    @NotNull
    private final ReportResources reportResources;

    public RNAFindingsChapter(@NotNull final OrangeRecord report, @NotNull final ReportResources reportResources) {
        this.report = report;
        this.reportResources = reportResources;
    }

    @NotNull
    @Override
    public String name() {
        return "RNA Findings";
    }

    @NotNull
    @Override
    public PageSize pageSize() {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph(name()).addStyle(reportResources.chapterTitleStyle()));

        addKeyQC(document);
        addExpressionTables(document);
        addRNAFusionTables(document);
        addNovelSpliceJunctionTables(document);
    }

    private void addKeyQC(@NotNull Document document) {
        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1 },
                new Cell[] { cells.createHeader("QC"), cells.createHeader("Total Fragments"), cells.createHeader("Non-Duplicate Fragments"),
                        cells.createHeader("Duplicate rate") });

        if (report.isofox() != null) {
            table.addCell(cells.createContent(report.isofox().summary().qcStatus()));
            table.addCell(cells.createContent(String.valueOf(report.isofox().summary().totalFragments())));

            long nonDuplicates = report.isofox().summary().totalFragments() - report.isofox().summary().duplicateFragments();
            table.addCell(cells.createContent(String.valueOf(nonDuplicates)));

            double duplicateRate = report.isofox().summary().duplicateFragments() / (double) report.isofox().summary().totalFragments();
            table.addCell(cells.createContent(PERCENTAGE_FORMAT.format(duplicateRate * 100)));
        } else {
            table.addCell(cells.createSpanningEntry(table, ReportResources.NOT_AVAILABLE));
        }

        document.add(new Tables(reportResources).createWrapping(table));
    }

    private void addExpressionTables(@NotNull Document document) {
        IsofoxRecord isofox = report.isofox();
        List<PurpleGeneCopyNumber> somaticGeneCopyNumbers = report.purple().allSomaticGeneCopyNumbers();

        List<GeneExpression> reportableHighExpression = isofox != null ? isofox.reportableHighExpression() : Lists.newArrayList();
        String titleHighExpression = "Genes with high expression (" + reportableHighExpression.size() + ")";
        document.add(ExpressionTable.build(titleHighExpression, contentWidth(), reportableHighExpression, false, somaticGeneCopyNumbers,
                reportResources));

        List<GeneExpression> reportableLowExpression = isofox != null ? isofox.reportableLowExpression() : Lists.newArrayList();
        String titleLowExpression = "Genes with low expression (" + reportableLowExpression.size() + ")";
        document.add(ExpressionTable.build(titleLowExpression, contentWidth(), reportableLowExpression, true, somaticGeneCopyNumbers,
                reportResources));
    }

    private void addRNAFusionTables(@NotNull Document document) {
        IsofoxRecord isofox = report.isofox();

        List<RnaFusion> reportableNovelKnownFusions = isofox != null ? isofox.reportableNovelKnownFusions() : Lists.newArrayList();
        String titleKnownFusions = "Known fusions detected in RNA and not in DNA (" + reportableNovelKnownFusions.size() + ")";
        document.add(RNAFusionTable.build(titleKnownFusions, contentWidth(), reportableNovelKnownFusions, reportResources));

        List<RnaFusion> reportableNovelPromiscuous = isofox != null ? isofox.reportableNovelPromiscuousFusions() : Lists.newArrayList();
        String titlePromiscuousFusions = "Promiscuous fusions detected in RNA and not in DNA (" + reportableNovelPromiscuous.size() + ")";
        document.add(RNAFusionTable.build(titlePromiscuousFusions, contentWidth(), reportableNovelPromiscuous, reportResources));
    }

    private void addNovelSpliceJunctionTables(@NotNull Document document) {
        IsofoxRecord isofox = report.isofox();

        List<NovelSpliceJunction> reportableSkippedExons = isofox != null ? isofox.reportableSkippedExons() : Lists.newArrayList();
        String titleSkippedExonJunctions =
                "Potentially interesting novel splice junctions - Skipped exons (" + reportableSkippedExons.size() + ")";
        document.add(NovelSpliceJunctionTable.build(titleSkippedExonJunctions, contentWidth(), reportableSkippedExons, reportResources));

        List<NovelSpliceJunction> reportableNovelExonsIntrons =
                isofox != null ? isofox.reportableNovelExonsIntrons() : Lists.newArrayList();
        String titleNovelExonIntronJunctions =
                "Potentially interesting novel splice junctions - Novel exon/intron (" + reportableNovelExonsIntrons.size() + ")";
        document.add(NovelSpliceJunctionTable.build(titleNovelExonIntronJunctions, contentWidth(), reportableNovelExonsIntrons,
                reportResources));
    }
}
