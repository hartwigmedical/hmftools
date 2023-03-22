package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.primitives.Doubles;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Expressions;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ExpressionTable {

    private static final Logger LOGGER = LogManager.getLogger(ExpressionTable.class);

    private ExpressionTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<GeneExpression> expressions, boolean sortAscending,
            @NotNull List<PurpleGeneCopyNumber> allSomaticGeneCopyNumbers) {
        if (expressions.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("Gene"), Cells.createHeader("Tumor CN"), Cells.createHeader("TPM"),
                        Cells.createHeader("Perc (Type)"), Cells.createHeader("FC (Type)"), Cells.createHeader("Perc (DB)"),
                        Cells.createHeader("FC (DB)") });

        // TODO Build the expression datamodel table prior to rendering.
        for (GeneExpression expression : sort(expressions, sortAscending)) {
            table.addCell(Cells.createContent(expression.geneName()));
            table.addCell(Cells.createContent(lookupTumorCN(allSomaticGeneCopyNumbers, expression.geneName())));
            table.addCell(Cells.createContent(Expressions.tpm(expression)));
            table.addCell(Cells.createContent(Expressions.percentileType(expression)));
            table.addCell(Cells.createContent(Expressions.foldChangeType(expression)));
            table.addCell(Cells.createContent(Expressions.percentileDatabase(expression)));
            table.addCell(Cells.createContent(Expressions.foldChangeDatabase(expression)));
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    private static String lookupTumorCN(@NotNull List<PurpleGeneCopyNumber> geneCopyNumbers, @NotNull String geneToFind) {
        PurpleGeneCopyNumber geneCopyNumber = findByGene(geneCopyNumbers, geneToFind);
        if (geneCopyNumber == null) {
            LOGGER.warn("Could not find gene copy number for '{}'", geneToFind);
            return ReportResources.NOT_AVAILABLE;
        }

        return String.valueOf(Math.round(Math.max(0, geneCopyNumber.minCopyNumber())));
    }

    @Nullable
    private static PurpleGeneCopyNumber findByGene(@NotNull List<PurpleGeneCopyNumber> geneCopyNumbers, @NotNull String geneToFind) {
        for (PurpleGeneCopyNumber geneCopyNumber : geneCopyNumbers) {
            if (geneCopyNumber.geneName().equals(geneToFind)) {
                return geneCopyNumber;
            }
        }

        return null;
    }

    @NotNull
    private static List<GeneExpression> sort(@NotNull List<GeneExpression> expressions, boolean sortDescending) {
        return expressions.stream().sorted((expression1, expression2) -> {
            if (sortDescending) {
                return Doubles.compare(expression1.percentileCohort(), expression2.percentileCohort());
            } else {
                return Doubles.compare(expression2.percentileCohort(), expression1.percentileCohort());
            }
        }).collect(Collectors.toList());
    }
}
