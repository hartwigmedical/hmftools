package com.hartwig.hmftools.orange.report.tables;

import static com.hartwig.hmftools.orange.report.tables.TableCommon.scaleTableWidths;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;

import org.junit.Test;

public class TableDataTest
{

    @Test
    public void testScalableTableData()
    {
        List<TableHeaderData> headers = Lists.newArrayList();

        List<List<String>> tableRowData = Lists.newArrayList();

        String dataLength10 = "0123456789";
        String dataLength5 = "01234";
        String dataLength15 = dataLength10 + dataLength5;
        String dataLength20 = dataLength10 + dataLength10;

        headers.add(new TableHeaderData(dataLength5, 0, false));
        headers.add(new TableHeaderData(dataLength5, 0, false));
        headers.add(new TableHeaderData(dataLength5, 0, false));

        tableRowData.add(List.of(dataLength10, dataLength5, dataLength10));
        tableRowData.add(List.of(dataLength5, dataLength10, dataLength10));
        tableRowData.add(List.of(dataLength10, dataLength10, dataLength10));

        // test 1: all 10, no padding
        int maxTableWidth = 30;
        scaleTableWidths(headers, tableRowData, maxTableWidth);

        assertEquals(10, headers.get(0).AdjustedWidth, 1);
        assertEquals(10, headers.get(1).AdjustedWidth, 1);
        assertEquals(10, headers.get(2).AdjustedWidth, 1);

        // test 2: adjust based on variable table data
        tableRowData.clear();

        tableRowData.add(List.of(dataLength15, dataLength5, dataLength10));
        tableRowData.add(List.of(dataLength5, dataLength5, dataLength5));
        tableRowData.add(List.of(dataLength10, dataLength5, dataLength5));

        // test 2: adjust to max data lengths
        scaleTableWidths(headers, tableRowData, maxTableWidth);

        assertEquals(15, headers.get(0).AdjustedWidth, 1);
        assertEquals(5, headers.get(1).AdjustedWidth, 1);
        assertEquals(10, headers.get(2).AdjustedWidth, 1);

        // test 3: adding buffer from excess
        scaleTableWidths(headers, tableRowData, maxTableWidth + 6);

        assertEquals(17, headers.get(0).AdjustedWidth, 1);
        assertEquals(7, headers.get(1).AdjustedWidth, 1);
        assertEquals(12, headers.get(2).AdjustedWidth, 1);

        // test 4: data columns wrapped
        headers.clear();

        headers.add(new TableHeaderData(dataLength5, 0, true));
        headers.add(new TableHeaderData(dataLength5, 0, false));
        headers.add(new TableHeaderData(dataLength5, 0, true));

        tableRowData.clear();

        tableRowData.add(List.of(dataLength20, dataLength5, dataLength10));
        tableRowData.add(List.of(dataLength15, dataLength5, dataLength20));
        tableRowData.add(List.of(dataLength10, dataLength20, dataLength5));

        maxTableWidth = 40;
        scaleTableWidths(headers, tableRowData, maxTableWidth);

        assertEquals(10, headers.get(0).AdjustedWidth, 1);
        assertEquals(20, headers.get(1).AdjustedWidth, 1);
        assertEquals(10, headers.get(2).AdjustedWidth, 1);

        // test 5: wrap headers while still keeping larger than row text
        headers.clear();

        headers.add(new TableHeaderData(dataLength5, 0, false));
        headers.add(new TableHeaderData(dataLength20, 0.5, false));
        headers.add(new TableHeaderData(dataLength5, 0, false));

        tableRowData.clear();

        tableRowData.add(List.of(dataLength5, dataLength5, dataLength20));
        tableRowData.add(List.of(dataLength5, dataLength5, dataLength5));
        tableRowData.add(List.of(dataLength20, dataLength5, dataLength5));

        maxTableWidth = 50;
        scaleTableWidths(headers, tableRowData, maxTableWidth);

        assertEquals(20, headers.get(0).AdjustedWidth, 1);
        assertEquals(10, headers.get(1).AdjustedWidth, 1);
        assertEquals(20, headers.get(2).AdjustedWidth, 1);
    }
}
