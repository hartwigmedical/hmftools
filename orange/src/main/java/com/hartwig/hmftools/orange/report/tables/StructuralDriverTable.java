package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.sv.linx.LinxDriver;
import com.hartwig.hmftools.orange.report.util.TableUtil;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class StructuralDriverTable {

    private StructuralDriverTable() {
    }

    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<LinxDriver> drivers) {
        if (drivers.isEmpty()) {
            return TableUtil.createEmptyTable(title, width);
        }

        Table table = TableUtil.createReportContentTable(width,
                new float[] { 1, 1, 4 },
                new Cell[] { TableUtil.createHeaderCell("Gene"), TableUtil.createHeaderCell("Event Type"),
                        TableUtil.createHeaderCell(Strings.EMPTY) });

        for (LinxDriver driver : sort(drivers)) {
            table.addCell(TableUtil.createContentCell(driver.gene()));
            table.addCell(TableUtil.createContentCell(driver.eventType()));
            table.addCell(TableUtil.createContentCell(Strings.EMPTY));
        }

        return TableUtil.createWrappingReportTable(table, title);
    }

    @NotNull
    private static List<LinxDriver> sort(@NotNull List<LinxDriver> drivers) {
        return drivers.stream().sorted((driver1, driver2) -> {
            if (driver1.gene().equals(driver2.gene())) {
                return driver1.eventType().compareTo(driver2.eventType());
            } else {
                return driver1.gene().compareTo(driver2.gene());
            }
        }).collect(Collectors.toList());
    }
}
