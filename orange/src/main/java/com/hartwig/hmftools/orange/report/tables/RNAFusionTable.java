package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.isofox.RnaFusion;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.jetbrains.annotations.NotNull;

public final class RNAFusionTable
{
    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<RnaFusion> fusions,
            @NotNull ReportResources reportResources)
    {
        if(fusions.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(width,
                new float[] { 2, 2, 2, 1, 2, 1, 2, 1 },
                new Cell[] { cells.createHeader("Name"), cells.createHeader("Pos (Up)"), cells.createHeader("Pos (Down)"),
                        cells.createHeader("SV Type"), cells.createHeader("Junction U/D"), cells.createHeader("Depth U/D"),
                        cells.createHeader("Frags (spl./re./disc.)"), cells.createHeader("Cohort freq") });

        for(RnaFusion fusion : sort(fusions))
        {
            table.addCell(cells.createContent(fusion.display()));
            table.addCell(cells.createContent(fusion.chromosomeUp() + ":" + fusion.positionUp()));
            table.addCell(cells.createContent(fusion.chromosomeDown() + ":" + fusion.positionDown()));
            table.addCell(cells.createContent(fusion.svType().toString()));
            table.addCell(cells.createContent(fusion.junctionTypeUp() + "/" + fusion.junctionTypeDown()));
            table.addCell(cells.createContent(fusion.depthUp() + "/" + fusion.depthDown()));
            table.addCell(cells.createContent(fusion.splitFragments() + "/" + fusion.realignedFrags() + "/" + fusion.discordantFrags()));
            table.addCell(cells.createContent(String.valueOf(fusion.cohortFrequency())));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    @NotNull
    private static List<RnaFusion> sort(@NotNull List<RnaFusion> fusions)
    {
        return fusions.stream().sorted((fusion1, fusion2) ->
        {
            String locationUp1 = Chromosomes.zeroPrefixed(fusion1.chromosomeUp());
            String locationUp2 = Chromosomes.zeroPrefixed(fusion2.chromosomeUp());

            if(locationUp1.equals(locationUp2))
            {
                return Integer.compare(fusion1.positionUp(), fusion2.positionUp());
            }
            else
            {
                return locationUp1.compareTo(locationUp2);
            }
        }).collect(Collectors.toList());
    }
}
