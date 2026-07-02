package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentileField;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Expressions;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.BaseTable;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatTpmField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_GENE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TPM;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.stringArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createEmptyTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.toPercentages;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import org.apache.logging.log4j.util.Strings;

public final class ExpressionTable
{
    private static String COL_PERCENTILE = "Percentile";
    private static String COL_FOLD_CHANGE = "Fold Change";

    public static BaseTable build(final DocumentContext docCtx,
            final String title, float width, final List<GeneExpression> expressions, boolean sortAscending,
            final ReportResources reportResources) throws IOException
    {
        if(expressions.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<String> headers = Lists.newArrayList();

        addEntry(widths, headers, 1, COL_GENE);
        addEntry(widths, headers, 1, COL_TPM);
        addEntry(widths, headers, 1, COL_PERCENTILE);
        addEntry(widths, headers, 1, COL_FOLD_CHANGE);
        addEntry(widths, headers, 3, Strings.EMPTY);

        BaseTable table = createStandardTable(docCtx, title, width, intToFloatArray(widths), stringArray(headers), reportResources);
        float[] pcts = toPercentages(intToFloatArray(widths));

        for(GeneExpression expression : sort(expressions, sortAscending))
        {
            List<String> rowValues = Lists.newArrayList();
            rowValues.add(expression.gene());
            rowValues.add(formatTpmField(expression.tpm()));

            if(expression.percentileCancer() != null && expression.medianTpmCancer() != null)
            {
                rowValues.add(formatPercentileField(expression.percentileCancer()));
                rowValues.add(Expressions.formatFoldChangeCancer(expression));
            }
            else
            {
                rowValues.add(formatPercentileField(expression.percentileCohort()));
                rowValues.add(Expressions.formatFoldChange(expression));
            }

            rowValues.add(Strings.EMPTY);
            cells.addRow(table, pcts, rowValues);
        }

        return table;
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
