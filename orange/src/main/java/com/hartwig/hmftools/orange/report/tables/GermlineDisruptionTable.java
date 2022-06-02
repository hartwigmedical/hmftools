package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.sv.linx.LinxGermlineSv;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class GermlineDisruptionTable {

    private GermlineDisruptionTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<LinxGermlineSv> disruptions) {
        if (disruptions.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("Gene"), Cells.createHeader("Pos (Start)"), Cells.createHeader("Pos (End)"),
                        Cells.createHeader("Type") });

        for (LinxGermlineSv disruption  : sort(disruptions)) {
            table.addCell(Cells.createContent(disruption.GeneName));
            table.addCell(Cells.createContent(disruption.ChromosomeStart + ":" + disruption.PositionStart));
            table.addCell(Cells.createContent(disruption.ChromosomeEnd + ":" + disruption.PositionEnd));
            table.addCell(Cells.createContent(disruption.Type.toString()));
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    private static List<LinxGermlineSv> sort(@NotNull List<LinxGermlineSv> disruptions) {
        return disruptions.stream().sorted((disruption1, disruption2) -> {
            int geneCompare = disruption1.GeneName.compareTo(disruption2.GeneName);
            if (geneCompare != 0) {
                return geneCompare;
            }

            int chrStartCompare = disruption1.ChromosomeStart.compareTo(disruption2.ChromosomeStart);
            if (chrStartCompare != 0) {
                return chrStartCompare;
            }

            return Integer.compare(disruption1.PositionStart, disruption2.PositionEnd);
        }).collect(Collectors.toList());
    }
}
