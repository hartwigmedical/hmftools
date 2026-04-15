package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentage;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.stringArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.toPercentages;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.isofox.RnaQCStatus;
import com.hartwig.hmftools.datamodel.isofox.RnaStatistics;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.BaseTable;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

public final class RnaStatisticsTable
{
    private static String COL_QC = "QC";
    private static String COL_TOTAL_FRAGS = "Total Fragments";
    private static String COL_DUP_RATE = "Duplicate Rate";
    private static String COL_SPLICED_RATE = "Spliced Rate";
    private static String COL_UNSPLICED_RATE = "Unspliced Rate";
    private static String COL_ALT_RATE = "Alt-sliced Rate";
    private static String COL_CHIMERIC_RATE = "Chimeric Rate";

    public static BaseTable build(final DocumentContext docCtx,
            final String title, float width, final RnaStatistics rnaStatistics, final ReportResources reportResources) throws IOException
    {
        Cells cells = new Cells(reportResources);

        List<Integer> widths = Lists.newArrayList();
        List<String> headers = Lists.newArrayList();

        addEntry(widths, headers, 1, COL_QC);
        addEntry(widths, headers, 1, COL_TOTAL_FRAGS);
        addEntry(widths, headers, 1, COL_DUP_RATE);
        addEntry(widths, headers, 1, COL_SPLICED_RATE);
        addEntry(widths, headers, 1, COL_UNSPLICED_RATE);
        addEntry(widths, headers, 1, COL_ALT_RATE);
        addEntry(widths, headers, 1, COL_CHIMERIC_RATE);

        BaseTable table = createStandardTable(docCtx, title, width, intToFloatArray(widths), stringArray(headers), reportResources);
        float[] pcts = toPercentages(intToFloatArray(widths));

        StringJoiner qcSj = new StringJoiner(", ");
        for(RnaQCStatus status : rnaStatistics.qcStatus())
        {
            qcSj.add(status.name());
        }

        List<String> rowValues = Lists.newArrayList();
        rowValues.add(qcSj.toString());
        rowValues.add(String.valueOf(rnaStatistics.totalFragments()));

        double duplicateRate = rnaStatistics.duplicateFragments() / (double) rnaStatistics.totalFragments();
        rowValues.add(formatPercentage(duplicateRate));

        rowValues.add(formatPercentage(rnaStatistics.splicedFragmentPerc()));
        rowValues.add(formatPercentage(rnaStatistics.unsplicedFragmentPerc()));
        rowValues.add(formatPercentage(rnaStatistics.altFragmentPerc()));
        rowValues.add(formatPercentage(rnaStatistics.chimericFragmentPerc()));

        cells.addRow(table, pcts, rowValues);

        return table;
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
}
