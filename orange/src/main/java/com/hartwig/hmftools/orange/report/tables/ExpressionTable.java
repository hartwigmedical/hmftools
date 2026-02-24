package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentileField;
import static com.hartwig.hmftools.orange.report.ReportResources.formatSingleDigitDecimal;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
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
import static com.hartwig.hmftools.orange.report.ReportResources.formatTpmField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_GENE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TPM;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.floatArray;

public final class ExpressionTable
{
    public static Table build(final String title, float width, final List<GeneExpression> expressions, boolean sortAscending,
            final List<PurpleGeneCopyNumber> allSomaticGeneCopyNumbers, final ReportResources reportResources)
    {
        if(expressions.isEmpty())
        {
            return new Tables(reportResources).createEmpty(title, width);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 1, COL_GENE);
        addEntry(cells, widths, cellEntries, 1, COL_TPM);
        addEntry(cells, widths, cellEntries, 1, "Percentile");
        addEntry(cells, widths, cellEntries, 1, "Fold Change");

        Table table = Tables.createContent(width, floatArray(widths), cellArray(cellEntries));
        
        for(GeneExpression expression : sort(expressions, sortAscending))
        {
            table.addCell(cells.createContent(expression.gene()));
            // table.addCell(cells.createContent(lookupTumorCN(allSomaticGeneCopyNumbers, expression.gene())));
            table.addCell(cells.createContent(formatTpmField(expression.tpm())));

            if(expression.percentileCancer() != null && expression.medianTpmCancer() != null)
            {
                table.addCell(cells.createContent(formatPercentileField(expression.percentileCancer())));
                table.addCell(cells.createContent(Expressions.formatFoldChangeCancer(expression)));

            }
            else
            {
                table.addCell(cells.createContent(formatPercentileField(expression.percentileCohort())));
                table.addCell(cells.createContent(Expressions.formatFoldChange(expression)));
            }
        }

        return new Tables(reportResources).createWrapping(table, title);
    }

    private static String lookupTumorCN(final List<PurpleGeneCopyNumber> geneCopyNumbers, final String geneToFind)
    {
        PurpleGeneCopyNumber geneCopyNumber = findByGene(geneCopyNumbers, geneToFind);
        if(geneCopyNumber == null)
        {
            LOGGER.warn("Could not find gene copy number for '{}'", geneToFind);
            return ReportResources.NOT_AVAILABLE;
        }

        return formatSingleDigitDecimal(Math.max(0, geneCopyNumber.minCopyNumber()));
    }

    private static PurpleGeneCopyNumber findByGene(final List<PurpleGeneCopyNumber> geneCopyNumbers, final String geneToFind)
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

    private static List<GeneExpression> sort(final List<GeneExpression> expressions, boolean sortDescending)
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
