package com.hartwig.hmftools.orange.report.tables;

import java.text.DecimalFormat;
import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.datamodel.BreakendEntry;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class BreakendTable {

    private static final DecimalFormat SINGLE_DIGIT = ReportResources.decimalFormat("#0.0");

    private BreakendTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<BreakendEntry> breakends) {
        if (breakends.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("Location"), Cells.createHeader("Gene"), Cells.createHeader("Range"),
                        Cells.createHeader("Type"), Cells.createHeader("Cluster ID"), Cells.createHeader("Junction CN"),
                        Cells.createHeader("Undisrupted CN") });

        for (BreakendEntry breakend : sort(breakends)) {
            table.addCell(Cells.createContent(breakend.location()));
            table.addCell(Cells.createContent(displayGene(breakend)));
            table.addCell(Cells.createContent(breakend.range()));
            table.addCell(Cells.createContent(breakend.type().toString()));
            table.addCell(Cells.createContent(String.valueOf(breakend.clusterId())));
            table.addCell(Cells.createContent(SINGLE_DIGIT.format(breakend.junctionCopyNumber())));
            table.addCell(Cells.createContent(SINGLE_DIGIT.format(breakend.undisruptedCopyNumber())));
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    public static List<BreakendEntry> sort(@NotNull List<BreakendEntry> breakends) {
        return breakends.stream().sorted((breakend1, breakend2) -> {
            String location1 = Chromosomes.zeroPrefixed(breakend1.location());
            String location2 = Chromosomes.zeroPrefixed(breakend2.location());

            int locationCompare = location1.compareTo(location2);
            if (locationCompare != 0) {
                return locationCompare;
            }

            int geneCompare = breakend1.gene().compareTo(breakend2.gene());
            if (geneCompare != 0) {
                return geneCompare;
            }

            return breakend1.exonUp() - breakend2.exonUp();
        }).collect(Collectors.toList());
    }

    @NotNull
    private static String displayGene(@NotNull BreakendEntry breakend) {
        String addon = Strings.EMPTY;
        if (!breakend.canonical()) {
            addon = " (alt)";
        }
        return breakend.gene() + addon;
    }
}
