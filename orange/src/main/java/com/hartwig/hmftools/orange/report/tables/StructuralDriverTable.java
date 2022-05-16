package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.sv.linx.LinxDriver;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
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
            return Tables.createEmpty(title, width);
        }

        Table table = Tables.createContent(width,
                new float[] { 1, 1, 4 },
                new Cell[] { Cells.createHeader("Gene"), Cells.createHeader("Event Type"), Cells.createHeader(Strings.EMPTY) });

        for (LinxDriver driver : sort(drivers)) {
            table.addCell(Cells.createContent(driver.gene()));
            table.addCell(Cells.createContent(driver.eventType()));
            table.addCell(Cells.createContent(Strings.EMPTY));
        }

        return Tables.createWrapping(table, title);
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
