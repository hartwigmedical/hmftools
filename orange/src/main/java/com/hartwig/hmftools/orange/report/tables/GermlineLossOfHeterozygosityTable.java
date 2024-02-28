package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.IsofoxRecord;
import com.hartwig.hmftools.datamodel.purple.GeneProportion;
import com.hartwig.hmftools.datamodel.purple.PurpleLossOfHeterozygosity;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Chromosomes;
import com.hartwig.hmftools.orange.report.interpretation.Expressions;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class GermlineLossOfHeterozygosityTable
{
    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<PurpleLossOfHeterozygosity> lossOfHeterozygosities,
            @Nullable IsofoxRecord isofox, @NotNull ReportResources reportResources)
    {
        if(lossOfHeterozygosities.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { cells.createHeader("Location"), cells.createHeader("Gene"),
                        cells.createHeader("Type"), cells.createHeader("CN"), cells.createHeader("TPM"), cells.createHeader("Perc (Type)"),
                        cells.createHeader("FC (Type)"), cells.createHeader("Perc (DB)"), cells.createHeader("FC (DB)") });

        for(PurpleLossOfHeterozygosity lossOfHeterozygosity : sort(lossOfHeterozygosities))
        {
            table.addCell(cells.createContent(lossOfHeterozygosity.chromosome() + lossOfHeterozygosity.chromosomeBand()));
            table.addCell(cells.createContent(displayGene(lossOfHeterozygosity)));
            table.addCell(cells.createContent(displayProportion(lossOfHeterozygosity.geneProportion())));
            table.addCell(cells.createContent(formatSingleDigitDecimal(lossOfHeterozygosity.minCopies())));

            GeneExpression expression = findExpressionForGene(isofox, lossOfHeterozygosity.gene());
            if(expression != null)
            {
                table.addCell(cells.createContent(Expressions.tpm(expression)));
                table.addCell(cells.createContent(Expressions.percentileType(expression)));
                table.addCell(cells.createContent(Expressions.foldChangeType(expression)));
                table.addCell(cells.createContent(Expressions.percentileDatabase(expression)));
                table.addCell(cells.createContent(Expressions.foldChangeDatabase(expression)));
            }
            else
            {
                table.addCell(cells.createContent(ReportResources.NOT_AVAILABLE));
                table.addCell(cells.createContent(ReportResources.NOT_AVAILABLE));
                table.addCell(cells.createContent(ReportResources.NOT_AVAILABLE));
                table.addCell(cells.createContent(ReportResources.NOT_AVAILABLE));
                table.addCell(cells.createContent(ReportResources.NOT_AVAILABLE));
            }
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    @Nullable
    private static GeneExpression findExpressionForGene(@Nullable IsofoxRecord isofox, @NotNull String geneToFind)
    {
        if(isofox == null)
        {
            return null;
        }

        return Expressions.findByGene(isofox.allGeneExpressions(), geneToFind);
    }

    @NotNull
    private static List<PurpleLossOfHeterozygosity> sort(@NotNull List<PurpleLossOfHeterozygosity> lossOfHeterozygosities)
    {
        return lossOfHeterozygosities.stream().sorted((loh1, loh2) ->
        {
            String location1 = Chromosomes.zeroPrefixed(loh1.chromosome() + loh1.chromosomeBand());
            String location2 = Chromosomes.zeroPrefixed(loh2.chromosome() + loh2.chromosomeBand());

            if(location1.equals(location2))
            {
                return loh1.gene().compareTo(loh2.gene());
            }
            else
            {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    private static String displayGene(@NotNull PurpleLossOfHeterozygosity deletion)
    {
        String addon = Strings.EMPTY;
        if(!deletion.isCanonical())
        {
            addon = " (alt)";
        }
        return deletion.gene() + addon;
    }

    @NotNull
    private static String displayProportion(@NotNull GeneProportion proportion)
    {
        return proportion.toString().toLowerCase().replaceAll("_", " ");
    }
}
