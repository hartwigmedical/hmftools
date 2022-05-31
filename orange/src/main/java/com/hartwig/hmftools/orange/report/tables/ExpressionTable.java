package com.hartwig.hmftools.orange.report.tables;

import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.rna.GeneExpression;
import com.hartwig.hmftools.orange.report.interpretation.Expressions;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ExpressionTable {

    private ExpressionTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<GeneExpression> expressions) {
        if (expressions.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("Chromosome"), Cells.createHeader("Region"), Cells.createHeader("Gene"),
                        Cells.createHeader("TPM"), Cells.createHeader("Perc (Type)"), Cells.createHeader("FC (Type)"),
                        Cells.createHeader("Perc (DB)"), Cells.createHeader("FC (DB)") });

        for (GeneExpression expression : sort(expressions)) {
            table.addCell(Cells.createContent(Strings.EMPTY));
            table.addCell(Cells.createContent(Strings.EMPTY));
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
    private static List<GeneExpression> sort(@NotNull List<GeneExpression> expressions) {
        return expressions.stream().sorted(Comparator.comparing(GeneExpression::geneName)).collect(Collectors.toList());
    }
}
