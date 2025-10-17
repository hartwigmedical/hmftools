package com.hartwig.hmftools.cobalt.ratio;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Before;
import org.junit.Test;

import tech.tablesaw.api.BooleanColumn;
import tech.tablesaw.api.DoubleColumn;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.Row;
import tech.tablesaw.api.StringColumn;
import tech.tablesaw.api.Table;

public class GcNormalisationTest
{
    private static final double EPSILON = 1e-5;

    @Before
    public void setUp()
    {
        // org.apache.logging.log4j.core.config.Configurator.setRootLevel(org.apache.logging.log4j.Level.TRACE);
    }

    @Test
    public void testGcNormaliser()
    {
        Table ratios = Table.create(
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION),
                DoubleColumn.create(CobaltColumns.READ_DEPTH),
                DoubleColumn.create(CobaltColumns.RATIO),
                DoubleColumn.create(CobaltColumns.READ_GC_CONTENT),
                BooleanColumn.create(CobaltColumns.IS_MAPPABLE),
                BooleanColumn.create(CobaltColumns.IS_AUTOSOME));

        addReadRatio(ratios, "chr1", 1001, 0.23, 0, 0.45, true);
        addReadRatio(ratios, "chr1", 2001,  50, 5, 0.451, true);
        addReadRatio(ratios, "chr1", 11001, 41, 4.0, 0.45, true);
        addReadRatio(ratios, "chr1", 12001, 180, 19, 0.501, true);
        addReadRatio(ratios, "chr2", 23001, 10, 1, 0.496, true);
        addReadRatio(ratios, "chr2", 24001, 20,2, 0.19, true); // gc bucket too low
        addReadRatio(ratios, "chr2", 25001, 30, 3, 0.70, true); // gc bucket too high
        addReadRatio(ratios, "chr3", 8001, 19,2, 0.45, false); // unmappable
        addReadRatio(ratios, "chrX", 7001, 21,2, 0.45, true); // allosome, not included in median calc

        (new ChromosomePositionCodec()).addEncodedChrPosColumn(ratios, false);

        ratios = new GcNormalizedRatioMapper(false).mapRatios(ratios);

        // System.out.println(ratios);

        assertEquals(6, ratios.rowCount());
        assertRatio(ratios, 0,1001, 0.0);
        assertRatio(ratios, 1,2001, 0.6896552);
        assertRatio(ratios, 2,11001, 0.5517241);
        assertRatio(ratios, 3,12001, 1.1793103);
        assertRatio(ratios, 4,23001, 0.062069);
        assertRatio(ratios, 5,7001, 0.275862);
    }

    private static void addReadRatio(Table table, String chromosome, int position, double depth, double ratio, double gcContent, boolean isMappable)
    {
        Row row = table.appendRow();
        row.setString(CobaltColumns.CHROMOSOME, chromosome);
        row.setInt(CobaltColumns.POSITION, position);
        row.setDouble(CobaltColumns.READ_DEPTH, depth);
        row.setDouble(CobaltColumns.RATIO, ratio);
        row.setDouble(CobaltColumns.READ_GC_CONTENT, gcContent);
        row.setBoolean(CobaltColumns.IS_MAPPABLE, isMappable);
        row.setBoolean(CobaltColumns.IS_AUTOSOME, HumanChromosome.fromString(chromosome).isAutosome());
    }

    private static void assertRatio(Table table, int rowIndex, int expectedPosition, double expectedRatio)
    {
        Row row = table.row(rowIndex);
        assertEquals(expectedPosition, row.getInt(CobaltColumns.POSITION));
        assertEquals(expectedRatio, row.getDouble(CobaltColumns.RATIO), EPSILON);
    }
}
