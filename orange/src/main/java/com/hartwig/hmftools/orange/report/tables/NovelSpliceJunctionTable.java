package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.isofox.NovelSpliceJunction;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class NovelSpliceJunctionTable
{
    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<NovelSpliceJunction> junctions,
            @NotNull ReportResources reportResources)
    {
        if(junctions.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(width,
                new float[] { 1, 1, 2, 2, 2, 2, 3, 1, 1 },
                new Cell[] { cells.createHeader("Gene"), cells.createHeader("Chr."), cells.createHeader("Junc (Start)"),
                        cells.createHeader("Junc (End)"), cells.createHeader("Type"), cells.createHeader("Depth S/E"),
                        cells.createHeader("Region S/E"), cells.createHeader("Frag Count"),
                        cells.createHeader("Cohort freq") });

        for(NovelSpliceJunction junction : sort(junctions))
        {
            table.addCell(cells.createContent(junction.gene()));
            table.addCell(cells.createContent(junction.chromosome()));
            table.addCell(cells.createContent(String.valueOf(junction.junctionStart())));
            table.addCell(cells.createContent(String.valueOf(junction.junctionEnd())));
            table.addCell(cells.createContent(junction.type().toString()));
            table.addCell(cells.createContent(junction.depthStart() + "/" + junction.depthEnd()));
            table.addCell(cells.createContent(junction.regionStart() + "/" + junction.regionEnd()));
            table.addCell(cells.createContent(String.valueOf(junction.fragmentCount())));
            table.addCell(cells.createContent(String.valueOf(junction.cohortFrequency())));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    @NotNull
    private static List<NovelSpliceJunction> sort(@NotNull List<NovelSpliceJunction> junctions)
    {
        return junctions.stream().sorted((junction1, junction2) ->
        {
            String locationUp1 = Chromosomes.zeroPrefixed(junction1.chromosome());
            String locationUp2 = Chromosomes.zeroPrefixed(junction2.chromosome());

            if(locationUp1.equals(locationUp2))
            {
                return Integer.compare(junction1.junctionStart(), junction2.junctionStart());
            }
            else
            {
                return locationUp1.compareTo(locationUp2);
            }
        }).collect(Collectors.toList());
    }
}
