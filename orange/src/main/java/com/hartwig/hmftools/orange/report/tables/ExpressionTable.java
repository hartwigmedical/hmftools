package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.primitives.Doubles;
import com.hartwig.hmftools.datamodel.purple.PurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Expressions;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ExpressionTable
{
    @NotNull
    public static Table build(@NotNull String title, float width, @NotNull List<GeneExpression> expressions, boolean sortAscending,
            @NotNull List<PurpleGeneCopyNumber> allSomaticGeneCopyNumbers, @NotNull ReportResources reportResources)
    {
        if(expressions.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);
        Table table = Tables.createContent(width,
                new float[] { 1, 1, 1, 1, 1, 1, 1 },
                new Cell[] { cells.createHeader("Gene"), cells.createHeader("Tumor CN"), cells.createHeader("TPM"),
                        cells.createHeader("Perc (Type)"), cells.createHeader("FC (Type)"), cells.createHeader("Perc (DB)"),
                        cells.createHeader("FC (DB)") });

        // TODO Build the expression datamodel table prior to rendering.
        for(GeneExpression expression : sort(expressions, sortAscending))
        {
            table.addCell(cells.createContent(expression.gene()));
            table.addCell(cells.createContent(lookupTumorCN(allSomaticGeneCopyNumbers, expression.gene())));
            table.addCell(cells.createContent(Expressions.tpm(expression)));
            table.addCell(cells.createContent(Expressions.percentileType(expression)));
            table.addCell(cells.createContent(Expressions.foldChangeType(expression)));
            table.addCell(cells.createContent(Expressions.percentileDatabase(expression)));
            table.addCell(cells.createContent(Expressions.foldChangeDatabase(expression)));
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    @NotNull
    private static String lookupTumorCN(@NotNull List<PurpleGeneCopyNumber> geneCopyNumbers, @NotNull String geneToFind)
    {
        PurpleGeneCopyNumber geneCopyNumber = findByGene(geneCopyNumbers, geneToFind);
        if(geneCopyNumber == null)
        {
            LOGGER.warn("Could not find gene copy number for '{}'", geneToFind);
            return ReportResources.NOT_AVAILABLE;
        }

        return formatSingleDigitDecimal(Math.max(0, geneCopyNumber.minCopyNumber()));
    }

    @Nullable
    private static PurpleGeneCopyNumber findByGene(@NotNull List<PurpleGeneCopyNumber> geneCopyNumbers, @NotNull String geneToFind)
    {
        for(PurpleGeneCopyNumber geneCopyNumber : geneCopyNumbers)
        {
            if(geneCopyNumber.gene().equals(geneToFind))
            {
                return geneCopyNumber;
            }
        }

        return null;
    }

    @NotNull
    private static List<GeneExpression> sort(@NotNull List<GeneExpression> expressions, boolean sortDescending)
    {
        return expressions.stream().sorted((expression1, expression2) ->
        {
            if(sortDescending)
            {
                return Doubles.compare(expression1.percentileCohort(), expression2.percentileCohort());
            }
            else
            {
                return Doubles.compare(expression2.percentileCohort(), expression1.percentileCohort());
            }
        }).collect(Collectors.toList());
    }
}
