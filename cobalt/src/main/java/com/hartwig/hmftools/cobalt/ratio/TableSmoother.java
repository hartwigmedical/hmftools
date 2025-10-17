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
        Table toSmooth = mInput.where(mInput.intColumn(0).isGreaterThan(0));
        boolean row0Presnt = toSmooth.rowCount() > mInput.rowCount();
        int rowCount = toSmooth.rowCount();
        if(rowCount < 3)
        {
            return toSmooth;
        }
        final String name0 = toSmooth.column(0).name();
        Table output = Table.create().addColumns(
                IntColumn.create(name0),
                DoubleColumn.create(toSmooth.column(1).name()),
                DoubleColumn.create(toSmooth.column(2).name())
        );
        toSmooth.sortOn(name0).rollingStream(3).forEach(rows -> addSmoothedRow(output, rows));

        if (row0Presnt)
        {
            Row newRow = output.appendRow();
            newRow.setInt(0, 0);
            newRow.setDouble(1, 1);
            newRow.setDouble(2, 1);
        }
        output.sortOn(name0);
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
