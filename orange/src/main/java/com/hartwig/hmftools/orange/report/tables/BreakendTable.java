package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.finding.BreakendEntry;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class BreakendTable
{
    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<BreakendEntry> breakends,
            @NotNull ReportResources reportResources)
    {
        if(breakends.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { cells.createHeader("Location"), cells.createHeader("Gene"), cells.createHeader("Range"),
                        cells.createHeader("Type"), cells.createHeader("Cluster ID"), cells.createHeader("Junction CN"),
                        cells.createHeader("Undisrupted CN") });

        for(BreakendEntry breakend : sort(breakends))
        {
            table.addCell(cells.createContent(breakend.location()));
            table.addCell(cells.createContent(displayGene(breakend)));
            table.addCell(cells.createContent(breakend.range()));
            table.addCell(cells.createContent(breakend.type().toString()));
            table.addCell(cells.createContent(String.valueOf(breakend.clusterId())));
            table.addCell(cells.createContent(formatSingleDigitDecimal(breakend.junctionCopyNumber())));
            table.addCell(cells.createContent(formatSingleDigitDecimal(breakend.undisruptedCopyNumber())));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    @NotNull
    private static List<BreakendEntry> sort(@NotNull List<BreakendEntry> breakends)
    {
        return breakends.stream().sorted((breakend1, breakend2) ->
        {
            String location1 = Chromosomes.zeroPrefixed(breakend1.location());
            String location2 = Chromosomes.zeroPrefixed(breakend2.location());

            int locationCompare = location1.compareTo(location2);
            if(locationCompare != 0)
            {
                return locationCompare;
            }

            int geneCompare = breakend1.gene().compareTo(breakend2.gene());
            if(geneCompare != 0)
            {
                return geneCompare;
            }

            return breakend1.exonUp() - breakend2.exonUp();
        }).collect(Collectors.toList());
    }

    @NotNull
    private static String displayGene(@NotNull BreakendEntry breakend)
    {
        String addon = Strings.EMPTY;
        if(!breakend.canonical())
        {
            addon = " (alt)";
        }
        return breakend.gene() + addon;
    }
}
