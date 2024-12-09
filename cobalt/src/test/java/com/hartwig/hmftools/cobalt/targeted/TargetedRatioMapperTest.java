package com.hartwig.hmftools.cobalt.targeted;

import static com.hartwig.hmftools.cobalt.CobaltTestUtils.assertDoubleEquals;

import static org.junit.Assert.assertEquals;

import java.util.Arrays;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.Chromosome;
import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.common.cobalt.ReadRatio;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import tech.tablesaw.api.*;

public class TargetedRatioMapperTest
{
    private static final Chromosome CHROMOSOME = new Chromosome("chr1", 10000);

    @Test
    public void testOnTargetRatio()
    {
        ChromosomePositionCodec chromosomePositionCodec = new ChromosomePositionCodec();

        final Table ratios = Table.create(
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION),
                DoubleColumn.create(CobaltColumns.RATIO),
                DoubleColumn.create(CobaltColumns.GC_CONTENT),
                IntColumn.create(CobaltColumns.GC_BUCKET),
                BooleanColumn.create(CobaltColumns.IS_MAPPABLE),
                BooleanColumn.create(CobaltColumns.IS_AUTOSOME));

        addReadRatio(ratios, 1001, 0, 45);
        addReadRatio(ratios, 2001, 0.5, 45);
        addReadRatio(ratios, 11001, 4.0, 45);
        addReadRatio(ratios, 12001, 19.5, 45);
        addReadRatio(ratios, 23001, 0, 45);

        chromosomePositionCodec.addEncodedChrPosColumn(ratios, false);

        Table targetEnrichmentRatios = Table.create(
                StringColumn.create(CobaltColumns.CHROMOSOME),
                IntColumn.create(CobaltColumns.POSITION),
                DoubleColumn.create(CobaltColumns.RELATIVE_ENRICHMENT),
                BooleanColumn.create("offTarget"));

        Row row = targetEnrichmentRatios.appendRow();
        row.setString(CobaltColumns.CHROMOSOME, CHROMOSOME.contig);
        row.setInt(CobaltColumns.POSITION, 2001);
        row.setDouble(CobaltColumns.RELATIVE_ENRICHMENT, 2.0);
        row.setBoolean("offTarget", false);

        row = targetEnrichmentRatios.appendRow();
        row.setString(CobaltColumns.CHROMOSOME, CHROMOSOME.contig);
        row.setInt(CobaltColumns.POSITION, 12001);
        row.setDouble(CobaltColumns.RELATIVE_ENRICHMENT, 10.0);
        row.setBoolean("offTarget", false);

        chromosomePositionCodec.addEncodedChrPosColumn(targetEnrichmentRatios, true);

        TargetedRatioMapper ratioMapper = new TargetedRatioMapper(targetEnrichmentRatios, chromosomePositionCodec);

        Table onTargetRatios = ratioMapper.onTargetRatios(ratios);

        assertEquals(2, onTargetRatios.rowCount());

        Row readRatio2001 = onTargetRatios.row(0);
        Row readRatio12001 = onTargetRatios.row(1);

        assertEquals(2001, readRatio2001.getInt(CobaltColumns.POSITION));

        // ratio = raw ratio / target enrichment / median of raw ratios that overlap with targeted

        // median of the unnormalized gc ratio is 10.0
        // so read ratio = 0.5 / 2.0 / 10 = 0.025
        assertDoubleEquals(0.025, readRatio2001.getDouble(CobaltColumns.RATIO));

        assertEquals(12001, readRatio12001.getInt(CobaltColumns.POSITION));

        // median of the unnormalized gc ratio is 10.0
        // so read ratio = 19.5 / 10.0 / 10 = 0.195
        assertDoubleEquals(0.195, readRatio12001.getDouble(CobaltColumns.RATIO));
    }

    @NotNull
    private static ListMultimap<Chromosome, ReadRatio> create(ReadRatio ... readRatios)
    {
        ListMultimap<Chromosome, ReadRatio> ratios = ArrayListMultimap.create();
        ratios.putAll(CHROMOSOME, Arrays.asList(readRatios));
        return ratios;
    }

    private static void addReadRatio(Table table, int position, double ratio, int gcBucket)
    {
        Row row = table.appendRow();
        row.setString(CobaltColumns.CHROMOSOME, CHROMOSOME.contig);
        row.setInt(CobaltColumns.POSITION, position);
        row.setDouble(CobaltColumns.RATIO, ratio);
        row.setDouble(CobaltColumns.GC_CONTENT, gcBucket / 100.0);
        row.setInt(CobaltColumns.GC_BUCKET, gcBucket);
        row.setBoolean(CobaltColumns.IS_MAPPABLE, true);
        row.setBoolean(CobaltColumns.IS_AUTOSOME, true);
    }
}
