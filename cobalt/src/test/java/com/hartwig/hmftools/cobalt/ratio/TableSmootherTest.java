package com.hartwig.hmftools.cobalt.ratio;

import static com.hartwig.hmftools.cobalt.CobaltTestUtils.assertDoubleEquals;

import static org.junit.Assert.assertEquals;

import org.junit.Before;
import org.junit.Test;

import tech.tablesaw.api.DoubleColumn;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.Row;
import tech.tablesaw.api.Table;

public class TableSmootherTest
{
    private static final String COL_0 = "gcBucket";
    private static final String COL_1 = "gcMedianCount";
    private static final String COL_2 = "windowCount";

    private Table input;
    private Table output;

    @Before
    public void setup()
    {
        input = Table.create("input",
                IntColumn.create(COL_0),
                DoubleColumn.create(COL_1),
                DoubleColumn.create(COL_2));
    }

    @Test
    public void rowsOutOfOrder()
    {
        addRow(4, 1.4, 14);
        addRow(3, 1.3, 13);
        addRow(2, 1.2, 12);
        addRow(1, 1.1, 11);

        smoothIt();
        assertEquals(2, output.rowCount());
        checkRow(0, 2, 1.2, 12.0);
        checkRow(1, 3, 1.3, 13.0);
    }

    @Test
    public void empty()
    {
        smoothIt();
    }

    @Test
    public void oneRow()
    {
        addRow(1, 1.0, 10.0);

        smoothIt();
        checkRow(0, 1, 1.0, 10.0);
    }

    @Test
    public void twoRows()
    {
        addRow(1, 1.0, 10);
        addRow(2, 1.2, 12);

        smoothIt();
        assertEquals(2, output.rowCount());
        checkRow(0, 1, 1.0, 10.0);
        checkRow(1, 2, 1.2, 12.0);
    }

    @Test
    public void threeRows()
    {
        addRow(1, 1.1, 11);
        addRow(2, 1.2, 12);
        addRow(3, 1.3, 13);

        smoothIt();
        assertEquals(1, output.rowCount());
        checkRow(0, 2, 1.2, 12.0);
    }

    @Test
    public void smoothedTest()
    {
        addRow(1, 1.0, 10);
        addRow(2, 1.3, 13);
        addRow(3, 1.0, 10);
        addRow(4, 1.3, 13);
        addRow(5, 1.6, 16);
        addRow(6, 1.0, 10);

        smoothIt();
        assertEquals(4, output.rowCount());
        checkRow(0, 2, 1.1, 11.0);
        checkRow(1, 3, 1.2, 12.0);
        checkRow(2, 4, 1.3, 13.0);
        checkRow(3, 5, 1.3, 13);
    }


    @Test
    public void gapTest()
    {
        addRow(1, 1.0, 10);
        addRow(2, 1.3, 13);
        addRow(4, 1.3, 13);
        addRow(5, 1.6, 16);
        addRow(6, 1.0, 10);

        smoothIt();
        assertEquals(3, output.rowCount());
        checkRow(0, 2, 1.2, 12.0);
        checkRow(1, 4, 1.4, 14.0);
        checkRow(2, 5, 1.3, 13.0);
    }

    private void checkRow(int index, int col0, double col1, double col2)
    {
        Row row = output.row(index);
        assertEquals(col0, row.getInt(COL_0));
        assertDoubleEquals(col1, row.getDouble(COL_1));
        assertDoubleEquals(col2, row.getDouble(COL_2));
    }

    private void smoothIt()
    {
        output = new TableSmoother(input).smoothed();
    }

    private void addRow(int col0, double col1, double col2)
    {
        Row row = input.appendRow();
        row.setInt(COL_0, col0);
        row.setDouble(COL_1, col1);
        row.setDouble(COL_2, col2);

    }
}
