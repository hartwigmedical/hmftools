package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentageOneDecimal;

import java.util.Comparator;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class MissedVariantLikelihoodTable
{
    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull Map<String, Double> significantGermlineMVLHPerGene,
            @NotNull ReportResources reportResources)
    {
        if(significantGermlineMVLHPerGene.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(width,
                new float[] { 2, 2, 1, 2, 2, 1, 2, 2, 1, 2, 2, 1 },
                new Cell[] { cells.createHeader("Gene"), cells.createHeader("MVLH"), cells.createHeader(Strings.EMPTY),
                        cells.createHeader("Gene"), cells.createHeader("MVLH"), cells.createHeader(Strings.EMPTY),
                        cells.createHeader("Gene"), cells.createHeader("MVLH"), cells.createHeader(Strings.EMPTY),
                        cells.createHeader("Gene"), cells.createHeader("MVLH"), cells.createHeader(Strings.EMPTY) });

        Set<String> genes = Sets.newTreeSet(Comparator.naturalOrder());
        genes.addAll(significantGermlineMVLHPerGene.keySet());
        for(String gene : genes)
        {
            double mvlh = significantGermlineMVLHPerGene.get(gene);
            table.addCell(cells.createContent(gene));
            table.addCell(cells.createContent(formatPercentageOneDecimal(mvlh)));
            table.addCell(cells.createContent(Strings.EMPTY));
        }

        // Make sure all rows are properly filled in case table is sparse.
        int entryCount = significantGermlineMVLHPerGene.size();
        if(entryCount % 4 != 0)
        {
            for(int i = 0; i < 12 - 3 * (entryCount % 4); i++)
            {
                table.addCell(cells.createContent(Strings.EMPTY));
            }
        }

        return new Tables(reportResources).createWrapping(table, title);
    }
}
