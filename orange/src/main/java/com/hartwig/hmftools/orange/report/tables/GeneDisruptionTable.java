package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.orange.algo.linx.GeneDisruption;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class GeneDisruptionTable {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#0.0");

    private GeneDisruptionTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<GeneDisruption> disruptions) {
        if (disruptions.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("Location"), Cells.createHeader("Gene"), Cells.createHeader("Range"),
                        Cells.createHeader("Type"), Cells.createHeader("Cluster ID"), Cells.createHeader("Junction CN"),
                        Cells.createHeader("Undisrupted CN") });

        for (GeneDisruption disruption : sort(disruptions)) {
            table.addCell(Cells.createContent(disruption.location()));
            table.addCell(Cells.createContent(displayGene(disruption)));
            table.addCell(Cells.createContent(disruption.range()));
            table.addCell(Cells.createContent(disruption.type()));
            table.addCell(Cells.createContent(String.valueOf(disruption.clusterId())));
            table.addCell(Cells.createContent(SINGLE_DIGIT.format(disruption.junctionCopyNumber())));
            table.addCell(Cells.createContent(SINGLE_DIGIT.format(disruption.undisruptedCopyNumber())));
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    public static List<GeneDisruption> sort(@NotNull List<GeneDisruption> disruptions) {
        return disruptions.stream().sorted((disruption1, disruption2) -> {
            String locationAndGene1 = Chromosomes.zeroPrefixed(disruption1.location()) + disruption1.gene();
            String locationAndGene2 = Chromosomes.zeroPrefixed(disruption2.location()) + disruption2.gene();

            if (locationAndGene1.equals(locationAndGene2)) {
                return disruption1.firstAffectedExon() - disruption2.firstAffectedExon();
            } else {
                return locationAndGene1.compareTo(locationAndGene2);
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    private static String displayGene(@NotNull GeneDisruption disruption) {
        String addon = Strings.EMPTY;
        if (!disruption.isCanonical()) {
            addon = " (alt)";
        }
        return disruption.gene() + addon;
    }
}
