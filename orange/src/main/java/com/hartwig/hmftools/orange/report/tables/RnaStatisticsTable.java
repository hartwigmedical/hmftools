package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentage;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.addEntry;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.cellArray;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.intToFloatArray;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.isofox.RnaQCStatus;
import com.hartwig.hmftools.datamodel.isofox.RnaStatistics;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;
import com.hartwig.hmftools.orange.report.util.Tables;
import com.itextpdf.layout.element.Cell;
import com.itextpdf.layout.element.Table;

public final class RnaStatisticsTable
{
    private static String COL_QC = "QC";
    private static String COL_TOTAL_FRAGS = "Total Fragments";
    private static String COL_DUP_RATE = "Duplicate Rate";
    private static String COL_SPLICED_RATE = "Spliced Rate";
    private static String COL_UNSPLICED_RATE = "Unspliced Rate";
    private static String COL_ALT_RATE = "Alt-sliced Rate";
    private static String COL_CHIMERIC_RATE = "Chimeric Rate";

    public static Table build(
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
}
