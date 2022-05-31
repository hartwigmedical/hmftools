package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class NovelSpliceJunctionTable {

    private NovelSpliceJunctionTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<NovelSpliceJunction> junctions) {
        if (junctions.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("Name"), Cells.createHeader("Pos (Up)"), Cells.createHeader("Pos (Down)"),
                        Cells.createHeader("SV Type"), Cells.createHeader("Junction Up/Down"), Cells.createHeader("Depth Up/Down"),
                        Cells.createHeader("Frag support (split/realigned/discordant)"), Cells.createHeader("Cohort freq") });

        for (NovelSpliceJunction junction : sort(junctions)) {
            // TODO
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    private static List<NovelSpliceJunction> sort(@NotNull List<NovelSpliceJunction> junctions) {
        return junctions.stream().sorted((junction1, junction2) -> {
            String locationUp1 = Chromosomes.zeroPrefixed(junction1.chromosome());
            String locationUp2 = Chromosomes.zeroPrefixed(junction2.chromosome());

            if (locationUp1.equals(locationUp2)) {
                return Integer.compare(junction1.junctionStart(), junction2.junctionStart());
            } else {
                return locationUp1.compareTo(locationUp2);
            }
        }).collect(Collectors.toList());
    }
}
