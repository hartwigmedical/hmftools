package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentage;
import static com.hartwig.hmftools.orange.report.ReportResources.formatPercentileField;

import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;
import com.hartwig.hmftools.datamodel.isofox.GeneExpression;
import com.hartwig.hmftools.datamodel.isofox.RnaStatistics;
import com.hartwig.hmftools.datamodel.isofox.RnaQCStatus;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.interpretation.Expressions;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

import static com.hartwig.hmftools.orange.report.ReportResources.formatTpmField;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_GENE;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.COL_TPM;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import org.apache.logging.log4j.util.Strings;

public final class ExpressionTable
{
    private static String COL_QC = "QC";
    private static String COL_TOTAL_FRAGS = "Total Fragments";
    private static String COL_DUP_RATE = "Duplicate Rate";
    private static String COL_SPLICED_RATE = "Spliced Rate";
    private static String COL_UNSPLICED_RATE = "Unspliced Rate";
    private static String COL_ALT_RATE = "Alt-sliced Rate";
    private static String COL_CHIMERIC_RATE = "Chimeric Rate";

    private static String COL_PERCENTILE = "Percentile";
    private static String COL_FOLD_CHANGE = "Fold Change";

    public static Table buildRnaSummary(
            final String title, float width, final RnaStatistics rnaStatistics, final ReportResources reportResources)
    {
        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<Cell> cellEntries = Lists.newArrayList();

        addEntry(cells, widths, cellEntries, 1, COL_QC);
        addEntry(cells, widths, cellEntries, 1, COL_TOTAL_FRAGS);
        addEntry(cells, widths, cellEntries, 1, COL_DUP_RATE);
        addEntry(cells, widths, cellEntries, 1, COL_SPLICED_RATE);
        addEntry(cells, widths, cellEntries, 1, COL_UNSPLICED_RATE);
        addEntry(cells, widths, cellEntries, 1, COL_ALT_RATE);
        addEntry(cells, widths, cellEntries, 1, COL_CHIMERIC_RATE);

        Table table = Tables.createContent(width, intToFloatArray(widths), cellArray(cellEntries));

        StringJoiner qcSj = new StringJoiner(", ");
        for(RnaQCStatus status : rnaStatistics.qcStatus())
        {
            qcSj.add(status.name());
        }

        table.addCell(cells.createContent(qcSj.toString()));
        table.addCell(cells.createContent(String.valueOf(rnaStatistics.totalFragments())));

        double duplicateRate = rnaStatistics.duplicateFragments() / (double) rnaStatistics.totalFragments();
        table.addCell(cells.createContent(formatPercentage(duplicateRate)));

        table.addCell(cells.createContent(formatPercentage(rnaStatistics.splicedFragmentPerc())));
        table.addCell(cells.createContent(formatPercentage(rnaStatistics.unsplicedFragmentPerc())));
        table.addCell(cells.createContent(formatPercentage(rnaStatistics.altFragmentPerc())));
        table.addCell(cells.createContent(formatPercentage(rnaStatistics.chimericFragmentPerc())));

        return new Tables(reportResources).createWrapping(table, title);
    }

    /*
    private void addQCWarningInCaseOfFail(Table table, Cells cells)
    {
        boolean isRnaFail = !mIsofoxRecord.summary().qcStatus().contains(RnaQCStatus.PASS);
        boolean isDnaFailNoTumor = PurpleQCInterpretation.isFailNoTumor(mPurpleRecord.fit().qc());

        if(isRnaFail || isDnaFailNoTumor)
        {
            String warning = isRnaFail ?
                    "The RNA QC status of this sample is not a pass. All presented RNA data should be interpreted with caution"
                    : "The DNA QC status of this sample is fail (no tumor). "
                            + "In addition to DNA findings, all RNA findings should be interpreted with caution";

            table.addCell(cells.createSpanningWarning(table, warning));
        }
    }
    */

    public static Table build(
            final String title, float width, final List<GeneExpression> expressions, boolean sortAscending, final ReportResources reportResources)
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
        addEntry(cells, widths, cellEntries, 1, COL_PERCENTILE);
        addEntry(cells, widths, cellEntries, 1, COL_FOLD_CHANGE);
        addEntry(cells, widths, cellEntries, 3, Strings.EMPTY);

        Table table = Tables.createContent(width, intToFloatArray(widths), cellArray(cellEntries));
        
        for(GeneExpression expression : sort(expressions, sortAscending))
        {
            table.addCell(cells.createContent(expression.gene()));
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

            table.addCell(cells.createContent(Strings.EMPTY));
        }

        return new Tables(reportResources).createWrapping(table, title);
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
