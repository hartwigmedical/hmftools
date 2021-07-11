package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.orange.report.util.TableUtil;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class GeneCopyNumberTable {

    private GeneCopyNumberTable() {
    }

    @Nullable
    public static Table build(@NotNull String title, @NotNull List<ReportableGainLoss> driverAmpsDels) {
        if (driverAmpsDels.isEmpty()) {
            return null;
        }

        Table table = TableUtil.createReportContentTable(new float[] { 1, 1, 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell("Chromosome"), TableUtil.createHeaderCell("Region"),
                        TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("Type"), TableUtil.createHeaderCell("Copies") });

        for (ReportableGainLoss gainLoss : sort(driverAmpsDels)) {
            table.addCell(TableUtil.createContentCell(gainLoss.chromosome()));
            table.addCell(TableUtil.createContentCell(gainLoss.chromosomeBand()));
            table.addCell(TableUtil.createContentCell(gainLoss.gene()));
            table.addCell(TableUtil.createContentCell(gainLoss.interpretation().display()));
            table.addCell(TableUtil.createContentCell(String.valueOf(gainLoss.copies())));
        }

        return TableUtil.createWrappingReportTable(table, title);
    }

    @NotNull
    private static List<ReportableGainLoss> sort(@NotNull List<ReportableGainLoss> reportableGainsAndLosses) {
        return reportableGainsAndLosses.stream().sorted((gainLoss1, gainLoss2) -> {
            String location1 = zeroPrefixed(gainLoss1.chromosome() + gainLoss1.chromosomeBand());
            String location2 = zeroPrefixed(gainLoss2.chromosome() + gainLoss2.chromosomeBand());

            if (location1.equals(location2)) {
                return gainLoss1.gene().compareTo(gainLoss2.gene());
            } else {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    private static String zeroPrefixed(@NotNull String location) {
        // First remove q or p arm if present.
        int armStart = location.indexOf("q");
        if (armStart < 0) {
            armStart = location.indexOf("p");
        }

        String chromosome = armStart > 0 ? location.substring(0, armStart) : location;

        try {
            int chromosomeIndex = Integer.parseInt(chromosome);
            if (chromosomeIndex < 10) {
                return "0" + location;
            } else {
                return location;
            }
        } catch (NumberFormatException exception) {
            // No need to prefix Y/X chromosomes
            return location;
        }
    }
}
