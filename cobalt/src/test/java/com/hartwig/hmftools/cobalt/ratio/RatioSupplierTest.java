package com.hartwig.hmftools.cobalt.ratio;

import static com.hartwig.hmftools.cobalt.CobaltColumns.ENCODED_CHROMOSOME_POS;
import static com.hartwig.hmftools.cobalt.CobaltColumns.TUMOR_GC_RATIO;
import static com.hartwig.hmftools.cobalt.CobaltColumns.TUMOR_READ_DEPTH;
import static com.hartwig.hmftools.cobalt.CobaltTestUtils.assertDoubleEquals;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

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

public class RatioSupplierTest
{

    private Table gcProfiles;
    private Table readDepths;
    private ChromosomePositionCodec chromosomePosCodec;
    private Table diploidRegions;

    @Before
    public void setup()
    {
        chromosomePosCodec = new ChromosomePositionCodec();
        // add some counts
        readDepths = Table.create("readDepths",
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION),
                DoubleColumn.create(CobaltColumns.READ_DEPTH),
                DoubleColumn.create(CobaltColumns.READ_GC_CONTENT));

        addReadDepth(readDepths, "chr1", 2001, 10.0);
        addReadDepth(readDepths, "chr2", 3001, 5.0);
        addReadDepth(readDepths, "chr2", 4001, 7.0);
        chromosomePosCodec.addEncodedChrPosColumn(readDepths, false);

        gcProfiles = Table.create("gcProfiles",
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION),
                DoubleColumn.create(CobaltColumns.GC_CONTENT),
                BooleanColumn.create(CobaltColumns.IS_MAPPABLE),
                BooleanColumn.create(CobaltColumns.IS_AUTOSOME));

        addGcProfile(gcProfiles, "chr1", 2001, 0.45, true);
        addGcProfile(gcProfiles, "chr2", 3001, 0.50, true);
        addGcProfile(gcProfiles, "chr2", 4001, 0.50, true);
        chromosomePosCodec.addEncodedChrPosColumn(gcProfiles, true);
    }

    @Test
    public void testTumorOnly() throws IOException
    {
        diploidRegions = Table.create("diploidRegions",
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
                gcProfiles, null, readDepths,
                chromosomePosCodec);

        Table ratios = ratioSupplier.tumorOnly(diploidRegions);

        assertEquals(3, ratios.rowCount());

        Row ratio = ratios.row(0);
        assertEquals(chromosomePosCodec.encodeChromosomePosition("chr1", 2001), ratio.getLong(ENCODED_CHROMOSOME_POS));
        assertDoubleEquals(10, ratio.getDouble(TUMOR_READ_DEPTH));
        assertDoubleEquals(0.95454545, ratio.getDouble(TUMOR_GC_RATIO));

        ratio = ratios.row(1);
        assertEquals(chromosomePosCodec.encodeChromosomePosition("chr2", 3001), ratio.getLong(ENCODED_CHROMOSOME_POS));
        assertDoubleEquals(5, ratio.getDouble(TUMOR_READ_DEPTH));

        ratio = ratios.row(2);
        assertEquals(chromosomePosCodec.encodeChromosomePosition("chr2", 4001), ratio.getLong(ENCODED_CHROMOSOME_POS));
        assertDoubleEquals(7, ratio.getDouble(TUMOR_READ_DEPTH));
        // tumorGCRatio must be -1 since this position is not in diploid bed file
        assertDoubleEquals(-1.0, ratio.getDouble(TUMOR_GC_RATIO));
    }

    @Test
    public void tumorOnlyNoDiploidData() throws IOException
    {
        diploidRegions = Table.create("diploidRegions",
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION));
        chromosomePosCodec.addEncodedChrPosColumn(diploidRegions, true);

        final RatioSupplier ratioSupplier = new RatioSupplier("TEST", "TEST", null,
                gcProfiles, null, readDepths,
                chromosomePosCodec);

        Table ratios = ratioSupplier.tumorOnly(diploidRegions);
        assertEquals(3, ratios.rowCount());
        Row ratio = ratios.row(2);
        // not masked out as no diploid regions file was supplied
        assertDoubleEquals(0.95454545, ratio.getDouble(TUMOR_GC_RATIO));
    }

    private static void addReadDepth(Table readCountTable, String chromosome, int position, double readDepth)
    {
        Row row = readCountTable.appendRow();
        row.setString(CobaltColumns.CHROMOSOME, chromosome);
        row.setInt(CobaltColumns.POSITION, position);
        row.setDouble(CobaltColumns.READ_DEPTH, readDepth);
        row.setDouble(CobaltColumns.READ_GC_CONTENT, 0.5);
    }

    @SuppressWarnings("SameParameterValue")
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
