package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class GeneCopyNumberTable {

    private GeneCopyNumberTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<ReportableGainLoss> gainLosses) {
        if (gainLosses.isEmpty()) {
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1 },
                new Cell[] { Cells.createHeader("Chromosome"), Cells.createHeader("Region"), Cells.createHeader("Gene"),
                        Cells.createHeader("Type"), Cells.createHeader("CN") });

        for (ReportableGainLoss gainLoss : sort(gainLosses)) {
            table.addCell(Cells.createContent(gainLoss.chromosome()));
            table.addCell(Cells.createContent(gainLoss.chromosomeBand()));
            table.addCell(Cells.createContent(gene(gainLoss)));
            table.addCell(Cells.createContent(gainLoss.interpretation().display()));
            table.addCell(Cells.createContent(String.valueOf(gainLoss.minCopies())));
        }

        return Tables.createWrapping(table, title);
    }

    @NotNull
    private static List<ReportableGainLoss> sort(@NotNull List<ReportableGainLoss> reportableGainsAndLosses) {
        return reportableGainsAndLosses.stream().sorted((gainLoss1, gainLoss2) -> {
            String location1 = ChromosomeUtil.zeroPrefixed(gainLoss1.chromosome() + gainLoss1.chromosomeBand());
            String location2 = ChromosomeUtil.zeroPrefixed(gainLoss2.chromosome() + gainLoss2.chromosomeBand());

            if (location1.equals(location2)) {
                return gainLoss1.gene().compareTo(gainLoss2.gene());
            } else {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    private static String gene(@NotNull ReportableGainLoss gainLoss) {
        String addon = Strings.EMPTY;
        if (!gainLoss.isCanonical()) {
            addon = " (alt)";
        }
        return gainLoss.gene() + addon;
    }
}
