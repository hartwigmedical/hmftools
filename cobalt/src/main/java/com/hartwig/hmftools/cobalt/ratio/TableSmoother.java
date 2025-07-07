package com.hartwig.hmftools.cobalt.ratio;

import java.util.Arrays;

import tech.tablesaw.api.DoubleColumn;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.Row;
import tech.tablesaw.api.Table;

public class TableSmoother
{
    private final Table mInput;

    public TableSmoother(final Table mInput)
    {
        this.mInput = mInput;
    }

    public Table smoothed()
    {
        int rowCount = mInput.rowCount();
        if(rowCount < 3)
        {
            return mInput;
        }
        final String name0 = mInput.column(0).name();
        Table output = Table.create().addColumns(
                IntColumn.create(name0),
                DoubleColumn.create(mInput.column(1).name()),
                DoubleColumn.create(mInput.column(2).name())
        );
        mInput.sortOn(name0).rollingStream(3).forEach(rows -> addSmoothedRow(output, rows));
        return output;
    }

    private void addSmoothedRow(Table output, Row[] rows)
    {
        int bucket = rows[1].getInt(0);
        double average1 = Arrays.stream(rows).mapToDouble(row -> row.getDouble(1)).average().orElseThrow();
        double average2 = Arrays.stream(rows).mapToDouble(row -> row.getDouble(2)).average().orElseThrow();
        Row newRow = output.appendRow();
        newRow.setInt(0, bucket);
        newRow.setDouble(1, average1);
        newRow.setDouble(2, average2);
    }
}
