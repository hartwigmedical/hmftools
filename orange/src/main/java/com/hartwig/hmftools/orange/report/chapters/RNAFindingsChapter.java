package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxInterpretedData;
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
    private final OrangeReport report;

    public RNAFindingsChapter(@NotNull final OrangeReport report) {
        this.report = report;
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
        document.add(new Paragraph(name()).addStyle(ReportResources.chapterTitleStyle()));

        addKeyQC(document);
        addExpressionTables(document);
        addRNAFusionTables(document);
        addNovelSpliceJunctionTables(document);
    }

    private void addKeyQC(@NotNull Document document) {
        Table table = Tables.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("QC"), Cells.createHeader("Total Fragments"), Cells.createHeader("Non-Duplicate Fragments"),
                        Cells.createHeader("Duplicate rate") });

        if (report.isofox() != null) {
            table.addCell(Cells.createContent(report.isofox().summary().qcStatus()));
            table.addCell(Cells.createContent(String.valueOf(report.isofox().summary().totalFragments())));

            long nonDuplicates = report.isofox().summary().totalFragments() - report.isofox().summary().duplicateFragments();
            table.addCell(Cells.createContent(String.valueOf(nonDuplicates)));

            double duplicateRate = report.isofox().summary().duplicateFragments() / (double) report.isofox().summary().totalFragments();
            table.addCell(Cells.createContent(PERCENTAGE_FORMAT.format(duplicateRate * 100)));
        } else {
            table.addCell(Cells.createSpanningEntry(table, ReportResources.NOT_AVAILABLE));
        }

        document.add(Tables.createWrapping(table));
    }

    private void addExpressionTables(@NotNull Document document) {
        IsofoxInterpretedData isofox = report.isofox();
        List<GeneCopyNumber> somaticGeneCopyNumbers = report.purple().allSomaticGeneCopyNumbers();

        List<GeneExpression> reportableHighExpression = isofox != null ? isofox.reportableHighExpression() : Lists.newArrayList();
        String titleHighExpression = "Genes with high expression (" + reportableHighExpression.size() + ")";
        document.add(ExpressionTable.build(titleHighExpression, contentWidth(), reportableHighExpression, false, somaticGeneCopyNumbers));

        List<GeneExpression> reportableLowExpression = isofox != null ? isofox.reportableLowExpression() : Lists.newArrayList();
        String titleLowExpression = "Genes with low expression (" + reportableLowExpression.size() + ")";
        document.add(ExpressionTable.build(titleLowExpression, contentWidth(), reportableLowExpression, true, somaticGeneCopyNumbers));
    }

    private void addRNAFusionTables(@NotNull Document document) {
        IsofoxInterpretedData isofox = report.isofox();

        List<RnaFusion> reportableNovelKnownFusions = isofox != null ? isofox.reportableNovelKnownFusions() : Lists.newArrayList();
        String titleKnownFusions = "Known fusions detected in RNA and not in DNA (" + reportableNovelKnownFusions.size() + ")";
        document.add(RNAFusionTable.build(titleKnownFusions, contentWidth(), reportableNovelKnownFusions));

        List<RnaFusion> reportableNovelPromiscuous = isofox != null ? isofox.reportableNovelPromiscuousFusions() : Lists.newArrayList();
        String titlePromiscuousFusions = "Promiscuous fusions detected in RNA and not in DNA (" + reportableNovelPromiscuous.size() + ")";
        document.add(RNAFusionTable.build(titlePromiscuousFusions, contentWidth(), reportableNovelPromiscuous));
    }

    private void addNovelSpliceJunctionTables(@NotNull Document document) {
        IsofoxInterpretedData isofox = report.isofox();

        List<NovelSpliceJunction> reportableSkippedExons = isofox != null ? isofox.reportableSkippedExons() : Lists.newArrayList();
        String titleSkippedExonJunctions =
                "Potentially interesting novel splice junctions - Skipped exons (" + reportableSkippedExons.size() + ")";
        document.add(NovelSpliceJunctionTable.build(titleSkippedExonJunctions, contentWidth(), reportableSkippedExons));

        List<NovelSpliceJunction> reportableNovelExonsIntrons =
                isofox != null ? isofox.reportableNovelExonsIntrons() : Lists.newArrayList();
        String titleNovelExonIntronJunctions =
                "Potentially interesting novel splice junctions - Novel exon/intron (" + reportableNovelExonsIntrons.size() + ")";
        document.add(NovelSpliceJunctionTable.build(titleNovelExonIntronJunctions, contentWidth(), reportableNovelExonsIntrons));
    }
}
