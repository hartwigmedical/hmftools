package com.hartwig.hmftools.cobalt.ratio;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Test;

import tech.tablesaw.api.BooleanColumn;
import tech.tablesaw.api.DoubleColumn;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.Row;
import tech.tablesaw.api.StringColumn;
import tech.tablesaw.api.Table;

public class RatioSupplierTest
{
    @Test
    public void testTumorOnly() throws IOException
    {
        var chromosomePosCodec = new ChromosomePositionCodec();

        // add some counts
        final Table readCounts = Table.create("readCounts",
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION),
                IntColumn.create(CobaltColumns.READ_COUNT));

        addReadCount(readCounts, "chr1", 2001, 100);
        addReadCount(readCounts, "chr2", 3001, 50);
        addReadCount(readCounts, "chr2", 4001, 70);

        chromosomePosCodec.addEncodedChrPosColumn(readCounts, false);

        // gc profiles
        Table gcProfiles = Table.create("gcProfiles",
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION),
                DoubleColumn.create(CobaltColumns.GC_CONTENT),
                BooleanColumn.create(CobaltColumns.IS_MAPPABLE),
                BooleanColumn.create(CobaltColumns.IS_AUTOSOME));

        addGcProfile(gcProfiles, "chr1", 2001, 0.45, true);
        addGcProfile(gcProfiles, "chr2", 3001, 0.50, true);
        addGcProfile(gcProfiles, "chr2", 4001, 0.50, true);

        chromosomePosCodec.addEncodedChrPosColumn(gcProfiles, true);

        // diploid regions
        final Table diploidRegions = Table.create("diploidRegions",
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION));

        Row row = diploidRegions.appendRow();
        row.setString(CobaltColumns.CHROMOSOME, "chr1");
        row.setInt(CobaltColumns.POSITION, 2001);

        row = diploidRegions.appendRow();
        row.setString(CobaltColumns.CHROMOSOME, "chr2");
        row.setInt(CobaltColumns.POSITION, 3001);
        chromosomePosCodec.addEncodedChrPosColumn(diploidRegions, true);

        final RatioSupplier ratioSupplier = new RatioSupplier("TEST", "TEST", null,
                gcProfiles, null, readCounts,
                chromosomePosCodec);

        Table ratios = ratioSupplier.tumorOnly(diploidRegions);

        // chr2:4001 is not in diploid region, therefore does not have ratio
        assertEquals(2, ratios.rowCount());

        Row ratio = ratios.row(0);
        assertEquals(chromosomePosCodec.encodeChromosomePosition("chr1", 2001), ratio.getLong(CobaltColumns.ENCODED_CHROMOSOME_POS));

        ratio = ratios.row(1);
        assertEquals(chromosomePosCodec.encodeChromosomePosition("chr2", 3001), ratio.getLong(CobaltColumns.ENCODED_CHROMOSOME_POS));
    }

    private static void addReadCount(Table readCountTable, String chromosome, int position, int readCount)
    {
        Row row = readCountTable.appendRow();
        row.setString(CobaltColumns.CHROMOSOME, chromosome);
        row.setInt(CobaltColumns.POSITION, position);
        row.setInt(CobaltColumns.READ_COUNT, readCount);
    }

    private static void addGcProfile(Table table, String chromosome, int position, double gcContent, boolean isMappable)
    {
        Row row = table.appendRow();
        row.setString(CobaltColumns.CHROMOSOME, chromosome);
        row.setInt(CobaltColumns.POSITION, position);
        row.setDouble(CobaltColumns.GC_CONTENT, gcContent);
        row.setBoolean(CobaltColumns.IS_MAPPABLE, isMappable);
        row.setBoolean(CobaltColumns.IS_AUTOSOME, HumanChromosome.fromString(chromosome).isAutosome());
    }
}
