package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class RNAFusionTable {

    private RNAFusionTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<RnaFusion> fusions) {
        if (fusions.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("Name"), Cells.createHeader("Pos (Up)"), Cells.createHeader("Pos (Down)"),
                        Cells.createHeader("SV Type"), Cells.createHeader("Junction U/D"), Cells.createHeader("Depth U/D"),
                        Cells.createHeader("Frags (spl./realig./disc.)"), Cells.createHeader("Cohort freq") });

        for (RnaFusion fusion : sort(fusions)) {
            table.addCell(Cells.createContent(fusion.name()));
            table.addCell(Cells.createContent(fusion.chromosomeUp() + ":" + fusion.positionUp()));
            table.addCell(Cells.createContent(fusion.chromosomeDown() + ":" + fusion.positionDown()));
            table.addCell(Cells.createContent(fusion.svType()));
            table.addCell(Cells.createContent(fusion.junctionTypeUp() + "/" + fusion.junctionTypeDown()));
            table.addCell(Cells.createContent(fusion.depthUp() + "/" + fusion.depthDown()));
            table.addCell(Cells.createContent(fusion.splitFragments() + "/" + fusion.realignedFrags() + "/" + fusion.discordantFrags()));
            table.addCell(Cells.createContent(String.valueOf(fusion.cohortFrequency())));
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    private static List<RnaFusion> sort(@NotNull List<RnaFusion> fusions) {
        return fusions.stream().sorted((fusion1, fusion2) -> {
            String locationUp1 = Chromosomes.zeroPrefixed(fusion1.chromosomeUp());
            String locationUp2 = Chromosomes.zeroPrefixed(fusion2.chromosomeUp());

            if (locationUp1.equals(locationUp2)) {
                return Integer.compare(fusion1.positionUp(), fusion2.positionUp());
            } else {
                return locationUp1.compareTo(locationUp2);
            }
        }).collect(Collectors.toList());
    }
}
