package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.primitives.Doubles;
import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.orange.report.interpretation.Expressions;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class ExpressionTable {

    private ExpressionTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<GeneExpression> expressions, boolean sortAscending) {
        if (expressions.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("Gene"), Cells.createHeader("TPM"), Cells.createHeader("Perc (Type)"),
                        Cells.createHeader("FC (Type)"), Cells.createHeader("Perc (DB)"), Cells.createHeader("FC (DB)") });

        for (GeneExpression expression : sort(expressions, sortAscending)) {
            table.addCell(Cells.createContent(expression.geneName()));
            table.addCell(Cells.createContent(Expressions.tpm(expression)));
            table.addCell(Cells.createContent(Expressions.percentileType(expression)));
            table.addCell(Cells.createContent(Expressions.foldChangeType(expression)));
            table.addCell(Cells.createContent(Expressions.percentileDatabase(expression)));
            table.addCell(Cells.createContent(Expressions.foldChangeDatabase(expression)));
        }

        return Tables.createWrapping(table, title);
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
