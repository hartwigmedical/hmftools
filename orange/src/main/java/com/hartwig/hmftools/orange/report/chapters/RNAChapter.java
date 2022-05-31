package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;

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

public class RNAChapter implements ReportChapter {

    private static final DecimalFormat PERCENTAGE_FORMAT = ReportResources.decimalFormat("#'%'");

    @NotNull
    private final OrangeReport report;

    public RNAChapter(@NotNull final OrangeReport report) {
        this.report = report;
    }

    @NotNull
    @Override
    public String name() {
        return "RNA";
    }

    @NotNull
    @Override
    public PageSize pageSize() {
        return PageSize.A4;
    }

    @Override
    public void render(@NotNull final Document document) {
        document.add(new Paragraph(name()).addStyle(ReportResources.chapterTitleStyle()));

        IsofoxInterpretedData isofox = report.isofox();
        if (isofox == null) {
            document.add(new Paragraph(ReportResources.NOT_AVAILABLE).addStyle(ReportResources.tableContentStyle()));
        } else {
            addKeyQC(document, isofox);
            addExpressionTables(document, isofox);
            addRNAFusionTables(document, isofox);
            addNovelSpliceJunctionTables(document, isofox);
        }
    }

    private void addKeyQC(@NotNull Document document, @NotNull IsofoxInterpretedData isofox) {
        Table table = Tables.createContent(contentWidth(),
                new float[] { 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("QC"), Cells.createHeader("Total Fragments"), Cells.createHeader("Non-Duplicate Fragments"),
                        Cells.createHeader("Duplicate rate") });

        table.addCell(Cells.createContent(isofox.summary().qcStatus()));
        table.addCell(Cells.createContent(String.valueOf(isofox.summary().totalFragments())));
        table.addCell(Cells.createContent(String.valueOf(isofox.summary().totalFragments() - isofox.summary().duplicateFragments())));

        double duplicateRate = isofox.summary().duplicateFragments() / (double) isofox.summary().totalFragments();
        table.addCell(Cells.createContent(PERCENTAGE_FORMAT.format(duplicateRate * 100)));

        document.add(Tables.createWrapping(table));
    }

    private void addExpressionTables(@NotNull Document document, @NotNull IsofoxInterpretedData isofox) {
        String titleHighExpression = "Genes with high expression (" + isofox.reportableHighExpression().size() + ")";
        document.add(ExpressionTable.build(titleHighExpression, contentWidth(), isofox.reportableHighExpression(), false));

        String titleLowExpression = "Genes with low expression (" + isofox.reportableLowExpression().size() + ")";
        document.add(ExpressionTable.build(titleLowExpression, contentWidth(), isofox.reportableLowExpression(), true));
    }

    private void addRNAFusionTables(@NotNull Document document, @NotNull IsofoxInterpretedData isofox) {
        String titleKnownFusions = "Known fusions detected in RNA and not in DNA (" + isofox.reportableNovelKnownFusions().size() + ")";
        document.add(RNAFusionTable.build(titleKnownFusions, contentWidth(), isofox.reportableNovelKnownFusions()));

        String titlePromiscuousFusions =
                "Promiscuous fusions detected in RNA and not in DNA (" + isofox.reportableNovelPromiscuousFusions().size() + ")";
        document.add(RNAFusionTable.build(titlePromiscuousFusions, contentWidth(), isofox.reportableNovelPromiscuousFusions()));
    }

    private void addNovelSpliceJunctionTables(@NotNull Document document, @NotNull IsofoxInterpretedData isofox) {
        String titleSkippedExonJunctions =
                "Potentially interesting novel splice junctions - Skipped exons (" + isofox.reportableSkippedExons().size() + ")";
        document.add(NovelSpliceJunctionTable.build(titleSkippedExonJunctions, contentWidth(), isofox.reportableSkippedExons()));

        String titleNovelExonIntronJunctions =
                "Potentially interesting novel splice junctions - Novel exon/intron (" + isofox.reportableNovelExonsIntrons().size() + ")";
        document.add(NovelSpliceJunctionTable.build(titleNovelExonIntronJunctions, contentWidth(), isofox.reportableNovelExonsIntrons()));
    }
}
