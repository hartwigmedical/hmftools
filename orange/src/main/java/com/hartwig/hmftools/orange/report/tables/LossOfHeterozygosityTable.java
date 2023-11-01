package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class LossOfHeterozygosityTable
{
    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<PurpleGeneCopyNumber> lohGenes,
            @NotNull ReportResources reportResources)
    {
        if(lohGenes.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 3 },
                new Cell[] { cells.createHeader("Location"), cells.createHeader("Gene"), cells.createHeader("Tumor MACN"),
                        cells.createHeader("Tumor CN"), cells.createHeader(Strings.EMPTY) });

        for(PurpleGeneCopyNumber lohGene : sort(lohGenes))
        {
            table.addCell(cells.createContent(lohGene.chromosome() + lohGene.chromosomeBand()));
            table.addCell(cells.createContent(lohGene.gene()));
            table.addCell(cells.createContent(formatSingleDigitDecimal(Math.max(0, lohGene.minMinorAlleleCopyNumber()))));
            table.addCell(cells.createContent(formatSingleDigitDecimal(Math.max(0, lohGene.minCopyNumber()))));
            table.addCell(cells.createContent(Strings.EMPTY));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    @NotNull
    private static List<PurpleGeneCopyNumber> sort(@NotNull List<PurpleGeneCopyNumber> lohGenes)
    {
        return lohGenes.stream().sorted((geneCopyNumber1, geneCopyNumber2) ->
        {
            String location1 = Chromosomes.zeroPrefixed(geneCopyNumber1.chromosome() + geneCopyNumber1.chromosomeBand());
            String location2 = Chromosomes.zeroPrefixed(geneCopyNumber2.chromosome() + geneCopyNumber2.chromosomeBand());

            if(location1.equals(location2))
            {
                return geneCopyNumber1.gene().compareTo(geneCopyNumber2.gene());
            }
            else
            {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }
}
