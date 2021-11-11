package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.orange.report.util.CellUtil;
import com.hartwig.hmftools.orange.report.util.TableUtil;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class GeneCopyNumberTable {

    private GeneCopyNumberTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<ReportableGainLoss> driverAmpsDels) {
        if (driverAmpsDels.isEmpty()) {
            return TableUtil.createEmpty(title, width);
        }

        Table table = TableUtil.createContent(width,
                new float[] { 1, 1, 1, 1, 1 },
                new Cell[] { CellUtil.createHeader("Chromosome"), CellUtil.createHeader("Region"), CellUtil.createHeader("Gene"),
                        CellUtil.createHeader("Type"), CellUtil.createHeader("CN") });

        for (ReportableGainLoss gainLoss : sort(driverAmpsDels)) {
            table.addCell(CellUtil.createContent(gainLoss.chromosome()));
            table.addCell(CellUtil.createContent(gainLoss.chromosomeBand()));
            table.addCell(CellUtil.createContent(gainLoss.gene()));
            table.addCell(CellUtil.createContent(gainLoss.interpretation().display()));
            table.addCell(CellUtil.createContent(String.valueOf(gainLoss.minCopies())));
        }

        return TableUtil.createWrapping(table, title);
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
}
