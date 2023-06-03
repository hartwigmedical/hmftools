package com.hartwig.hmftools.cobalt.ratio;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.cobalt.Chromosome;
import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;

import org.junit.Test;

import tech.tablesaw.api.BooleanColumn;
import tech.tablesaw.api.DoubleColumn;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.Row;
import tech.tablesaw.api.StringColumn;
import tech.tablesaw.api.Table;

public class GcNormalisationTest
{
    private static final Chromosome CHROMOSOME = new Chromosome("chr1", 10000);
    private static final double EPSILON = 1e-5;

    @Test
    public void testGcNormaliser()
    {
        Table ratios = Table.create(
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION),
                IntColumn.create(CobaltColumns.READ_COUNT),
                DoubleColumn.create(CobaltColumns.RATIO),
                IntColumn.create(CobaltColumns.GC_BUCKET),
                BooleanColumn.create(CobaltColumns.IS_MAPPABLE),
                BooleanColumn.create(CobaltColumns.IS_AUTOSOME));

        addReadRatio(ratios, 1001, 0, 45);
        addReadRatio(ratios, 2001, 5, 45);
        addReadRatio(ratios, 11001, 4.0, 45);
        addReadRatio(ratios, 12001, 19, 50);
        addReadRatio(ratios, 23001, 1, 50);

        (new ChromosomePositionCodec()).addEncodedChrPosColumn(ratios, false);

        var ratioBuilder = new GcNormalizedRatioBuilder(ratios, false);

        ratios = ratioBuilder.ratios();

        assertEquals(5, ratios.rowCount());
        assertRatio(ratios, 0,1001, 0.0);
        assertRatio(ratios, 1,2001, 0.86207);
        assertRatio(ratios, 3,12001, 1.31034);
    }

    private static void addReadRatio(Table table, int position, double ratio, int gcBucket)
    {
        Row row = table.appendRow();
        row.setString(CobaltColumns.CHROMOSOME, CHROMOSOME.contig);
        row.setInt(CobaltColumns.POSITION, position);
        row.setInt(CobaltColumns.READ_COUNT, 10);
        row.setDouble(CobaltColumns.RATIO, ratio);
        row.setInt(CobaltColumns.GC_BUCKET, gcBucket);
        row.setBoolean(CobaltColumns.IS_MAPPABLE, true);
        row.setBoolean(CobaltColumns.IS_AUTOSOME, true);
    }

    private static void assertRatio(Table table, int rowIndex, int expectedPosition, double expectedRatio)
    {
        Row row = table.row(rowIndex);
        assertEquals(expectedPosition, row.getInt(CobaltColumns.POSITION));
        assertEquals(expectedRatio, row.getDouble(CobaltColumns.RATIO), EPSILON);
    }
}
