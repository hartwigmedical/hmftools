package com.hartwig.hmftools.orange.report.tables;

import com.hartwig.hmftools.datamodel.purple.ChromosomalRearrangements;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ChromosomalRearrangementsTable
{
    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull ChromosomalRearrangements chromosomalRearrangements,
            @NotNull ReportResources reportResources)
    {
        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(width,
                new float[] { 2, 1, 3 },
                new Cell[] { cells.createHeader("Chromosomal rearrangement"), cells.createHeader("Detected?"),
                        cells.createHeader(Strings.EMPTY) });

        table.addCell(cells.createContent("1q trisomy"));
        table.addCell(cells.createContent(chromosomalRearrangements.hasTrisomy1q() ? "Yes" : "No"));
        table.addCell(cells.createContent(Strings.EMPTY));

        table.addCell(cells.createContent("1p19q co-deletion"));
        table.addCell(cells.createContent(chromosomalRearrangements.hasCodeletion1p19q() ? "Yes" : "No"));
        table.addCell(cells.createContent(Strings.EMPTY));

        return new Tables(reportResources).createWrapping(table, title);
    }
}
