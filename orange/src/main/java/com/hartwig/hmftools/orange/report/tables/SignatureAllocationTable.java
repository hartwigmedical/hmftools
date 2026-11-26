package com.hartwig.hmftools.orange.report.tables;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.sigs.SigResiduals.SIG_MISALLOCATED;
import static com.hartwig.hmftools.orange.report.ReportResources.MAX_TABLE_WIDTH;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.convertWidths;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.formatPercentage;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.scaleTableWidths;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createStandardTable;
import static com.hartwig.hmftools.orange.report.tables.TableCommon.createEmptyTable;

import java.util.List;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.sigs.SignatureAllocation;
import com.hartwig.hmftools.orange.report.ReportResources;
import com.hartwig.hmftools.orange.report.util.Cells;

import be.quodlibet.boxable.BaseTable;

import com.hartwig.hmftools.orange.report.DocumentContext;

import java.io.IOException;

public final class SignatureAllocationTable
{
    private static final String COL_SIGNATURE = "Signature";
    private static final String COL_ETIOLOGY = "Etiology";
    private static final String COL_ALLOCATION = "Allocation";
    private static final String COL_PERCENT = "Percent";

    public static BaseTable build(
            final DocumentContext docCtx, final List<SignatureAllocation> signatureAllocations, final String title, float width,
            final ReportResources reportResources)
            throws IOException
    {
        if(signatureAllocations.isEmpty())
        {
            return createEmptyTable(docCtx, title, width, reportResources);
        }

        List<TableHeaderData> headers = Lists.newArrayList();

        headers.add(new TableHeaderData(COL_SIGNATURE, 0, false));
        headers.add(new TableHeaderData(COL_ETIOLOGY, 0, false));
        headers.add(new TableHeaderData(COL_ALLOCATION, 0, false));
        headers.add(new TableHeaderData(COL_PERCENT, 0, false));

        List<List<String>> tableRowData = Lists.newArrayList();

        for(SignatureAllocation sigAllocation : sort(signatureAllocations))
        {
            List<String> rowData = Lists.newArrayList();

            tableRowData.add(rowData);

            if(sigAllocation.signature().equals(SIG_MISALLOCATED))
            {
                rowData.add("Mis-allocated");
                rowData.add("");
                rowData.add(format("%.0f", sigAllocation.allocation()));
                rowData.add("");
            }
            else
            {
                rowData.add(sigAllocation.signature());
                rowData.add(sigAllocation.etiology());
                rowData.add(format("%.0f", sigAllocation.allocation()));
                rowData.add(formatPercentage(sigAllocation.percent()));
            }
        }

        int tableWidth = MAX_TABLE_WIDTH - 100;
        scaleTableWidths(headers, tableRowData, tableWidth, reportResources);

        int columnCount = headers.size();
        String[] columnTexts = new String[columnCount];
        float[] columnWidths = new float[columnCount];
        float[] columnPercentages = new float[columnCount];

        convertWidths(headers, columnTexts, columnWidths, columnPercentages);

        BaseTable table = createStandardTable(docCtx, title, tableWidth, columnWidths, columnTexts, reportResources);

        Cells cells = new Cells(reportResources);

        for(int r = 0; r < tableRowData.size(); ++r)
        {
            List<String> rowData = tableRowData.get(r);
            cells.addRow(table, columnPercentages, rowData);
        }

        return table;
    }

    @VisibleForTesting
    static List<SignatureAllocation> sort(final List<SignatureAllocation> signatureAllocations)
    {
        return signatureAllocations.stream().sorted((allocation1, allocation2) ->
        {
            String sig1 = allocation1.signature().toLowerCase();
            String sig2 = allocation2.signature().toLowerCase();

            if(sig1.equalsIgnoreCase(SIG_MISALLOCATED))
            {
                return 1;
            }
            else if(sig2.equalsIgnoreCase(SIG_MISALLOCATED))
            {
                return -1;
            }


            return Double.compare(allocation2.allocation(), allocation1.allocation());
        }).collect(Collectors.toList());
    }
}
