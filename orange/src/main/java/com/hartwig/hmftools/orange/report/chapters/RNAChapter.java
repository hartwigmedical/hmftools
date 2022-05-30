package com.hartwig.hmftools.orange.report.chapters;

import java.text.DecimalFormat;

import com.hartwig.hmftools.common.isofox.IsofoxData;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.report.ReportResources;
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

        IsofoxData isofox = report.isofox();
        if (isofox == null) {
            document.add(new Paragraph(ReportResources.NOT_AVAILABLE).addStyle(ReportResources.tableContentStyle()));
        } else {
            addKeyQC(document, isofox);
            addHighExpressionGenes(document, isofox);
            addLowExpressionGenes(document, isofox);
            addNovelKnownRNAFusions(document, isofox);
            addNovelPromiscuousRNAFusions(document, isofox);
            addSkippedExonNovelSpliceJunctions(document, isofox);
            addNovelExonIntronNovelSpliceJunctions(document, isofox);
        }
    }

    private void addKeyQC(@NotNull Document document, @NotNull IsofoxData isofox) {
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

    private static void addHighExpressionGenes(@NotNull Document document, @NotNull IsofoxData isofox) {

    }

    private static void addLowExpressionGenes(@NotNull Document document, @NotNull IsofoxData isofox) {

    }

    private static void addNovelKnownRNAFusions(@NotNull Document document, @NotNull IsofoxData isofox) {

    }

    private static void addNovelPromiscuousRNAFusions(@NotNull Document document, @NotNull IsofoxData isofox) {

    }

    private static void addSkippedExonNovelSpliceJunctions(@NotNull Document document, @NotNull IsofoxData isofox) {

    }

    private static void addNovelExonIntronNovelSpliceJunctions(@NotNull Document document, @NotNull IsofoxData isofox) {

    }
}
