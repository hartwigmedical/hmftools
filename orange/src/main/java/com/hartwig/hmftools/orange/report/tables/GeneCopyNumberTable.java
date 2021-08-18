package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.copynumber.ReportableGainLoss;
import com.hartwig.hmftools.orange.report.util.GeneUtil;
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
            return TableUtil.createEmptyTable(title, width);
        }

        Table table = TableUtil.createReportContentTable(width, new float[] { 1, 1, 1, 1, 1 },
                new Cell[] { TableUtil.createHeaderCell("Chromosome"), TableUtil.createHeaderCell("Region"),
                        TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("Type"), TableUtil.createHeaderCell("CN") });

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
            String location1 = GeneUtil.zeroPrefixed(gainLoss1.chromosome() + gainLoss1.chromosomeBand());
            String location2 = GeneUtil.zeroPrefixed(gainLoss2.chromosome() + gainLoss2.chromosomeBand());

            if (location1.equals(location2)) {
                return gainLoss1.gene().compareTo(gainLoss2.gene());
            } else {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }
}
