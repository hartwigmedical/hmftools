package com.hartwig.hmftools.orange.report.tables;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.linx.LinxHomozygousDisruption;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class HomozygousDisruptionTable
{
    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<LinxHomozygousDisruption> homozygousDisruptions,
            @NotNull ReportResources reportResources)
    {
        if(homozygousDisruptions.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(width,
                new float[] { 1, 1, 4 },
                new Cell[] { cells.createHeader("Location"), cells.createHeader("Gene"), cells.createHeader(Strings.EMPTY) });

        for(LinxHomozygousDisruption homozygousDisruption : sort(homozygousDisruptions))
        {
            table.addCell(cells.createContent(homozygousDisruption.chromosome() + homozygousDisruption.chromosomeBand()));
            table.addCell(cells.createContent(gene(homozygousDisruption)));
            table.addCell(cells.createContent(Strings.EMPTY));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    @NotNull
    private static List<LinxHomozygousDisruption> sort(@NotNull List<LinxHomozygousDisruption> homozygousDisruptions)
    {
        return homozygousDisruptions.stream().sorted((disruption1, disruption2) ->
        {
            String location1 = Chromosomes.zeroPrefixed(disruption1.chromosome() + disruption1.chromosomeBand());
            String location2 = Chromosomes.zeroPrefixed(disruption2.chromosome() + disruption2.chromosomeBand());

            if(location1.equals(location2))
            {
                return disruption1.gene().compareTo(disruption2.gene());
            }
            else
            {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    private static String gene(@NotNull LinxHomozygousDisruption homozygousDisruption)
    {
        String addon = Strings.EMPTY;
        if(!homozygousDisruption.isCanonical())
        {
            addon = " (alt)";
        }
        return homozygousDisruption.gene() + addon;
    }
}
