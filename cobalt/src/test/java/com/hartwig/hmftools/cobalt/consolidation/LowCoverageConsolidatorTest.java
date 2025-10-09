package com.hartwig.hmftools.cobalt.consolidation;

import static com.hartwig.hmftools.cobalt.CobaltColumns.ENCODED_CHROMOSOME_POS;
import static com.hartwig.hmftools.cobalt.CobaltColumns.POSITION;
import static com.hartwig.hmftools.cobalt.CobaltColumns.RATIO;
import static com.hartwig.hmftools.cobalt.CobaltColumns.READ_GC_CONTENT;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;
import com.hartwig.hmftools.cobalt.calculations.BamRatio;
import com.hartwig.hmftools.cobalt.calculations.CalculationsTestBase;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Before;
import org.junit.Test;

import tech.tablesaw.api.BooleanColumn;
import tech.tablesaw.api.DoubleColumn;
import tech.tablesaw.api.IntColumn;
import tech.tablesaw.api.LongColumn;
import tech.tablesaw.api.Row;
import tech.tablesaw.api.StringColumn;
import tech.tablesaw.api.Table;

public class LowCoverageConsolidatorTest extends CalculationsTestBase
{
    private ChromosomePositionCodec Codec;

    @Before
    public void setup()
    {
        Codec = new ChromosomePositionCodec();
    }

    @Test
    public void ratioAndGcAreAveragedInConsolidation()
    {
        ListMultimap<Chromosome, BamRatio> ratios = ArrayListMultimap.create();
        for(int i = 0; i < 40; i++)
        {
            int position = i * 1000 + 1;
            double depth = 10 + i * 0.01;
            double gc = 0.40 + i * 0.01;
            ratios.put(_1, br(_1, position, depth, gc, true));
            ratios.put(_2, br(_2, position, depth, gc, true));
        }
        int consolidationCount = ResultsConsolidator.calcConsolidationCount(8.0);
        LowCoverageConsolidator consolidator = new LowCoverageConsolidator(consolidationCount);
        ListMultimap<Chromosome, BamRatio> result = consolidator.consolidate(ratios);

        Table asTable = toTable(ratios);

        Multimap<String, LowCovBucket> consolidated = LowCoverageRatioMapper.calcConsolidateBuckets(asTable, 8.0);
        assert consolidated != null;
        Table consolidatedTable = new LowCoverageRatioMapper(consolidated, Codec).mapRatios(asTable);
        ListMultimap<Chromosome, BamRatio> roundTripped = toRatios(consolidatedTable);

        assertEquals(2, roundTripped.keySet().size());
        List<BamRatio> ratios1 = roundTripped.get(_1);
        List<BamRatio> ratios2 = roundTripped.get(_2);
        assertEquals(4, ratios1.size());
        assertEquals(4, ratios2.size());
        assertEquals(4001, ratios1.get(0).Position);
        assertEquals(10.045, ratios1.get(0).ratio(), 0.0001);
        assertEquals(0.445, ratios1.get(0).gcContent(), 0.0001);
        assertEquals(4001, ratios2.get(0).Position);
        assertEquals(10.045, ratios2.get(0).ratio(), 0.0001);
        assertEquals(0.445, ratios2.get(0).gcContent(), 0.0001);
        assertEquals(14001, ratios1.get(1).Position);
        assertEquals(10.145, ratios1.get(1).ratio(), 0.0001);
        assertEquals(0.545, ratios1.get(1).gcContent(), 0.0001);
        assertEquals(24001, ratios1.get(2).Position);
        assertEquals(10.245, ratios1.get(2).ratio(), 0.0001);
        assertEquals(0.645, ratios1.get(2).gcContent(), 0.0001);
        assertEquals(35001, ratios1.get(3).Position);
        assertEquals(10.345, ratios1.get(3).ratio(), 0.0001);
        assertEquals(0.745, ratios1.get(3).gcContent(), 0.0001);

        assertEquals(roundTripped, result);
    }

    @Test
    public void maskedRatiosAreSkippedInConsolidation()
    {
        ListMultimap<Chromosome, BamRatio> ratios = ArrayListMultimap.create();
        for (int i = 0; i < 200; i++)
        {
            int position = i * 1000 + 1;
            if (i % 5 == 0)
            {
                ratios.put(_1, br(_1, position, 10.0, 0.5, true));
            }
            else
            {
                ratios.put(_1, br(_1, position, -1.0, 0.4, true));
            }
        }
        int consolidationCount = ResultsConsolidator.calcConsolidationCount(8.0);
        LowCoverageConsolidator consolidator = new LowCoverageConsolidator(consolidationCount);
        ListMultimap<Chromosome, BamRatio> result = consolidator.consolidate(ratios);

        Table asTable = toTable(ratios);

        Multimap<String, LowCovBucket> consolidated = LowCoverageRatioMapper.calcConsolidateBuckets(asTable, 8.0);
        assert consolidated != null;
        Table consolidatedTable = new LowCoverageRatioMapper(consolidated, Codec).mapRatios(asTable);
        ListMultimap<Chromosome, BamRatio> roundTripped = toRatios(consolidatedTable);

        assertEquals(1, roundTripped.keySet().size());
        List<BamRatio> ratios1 = roundTripped.get(_1);
        assertEquals(4, ratios1.size());
        assertEquals(23001, ratios1.get(0).Position);
        assertEquals(10.0, ratios1.get(0).ratio(), 0.0001);
        assertEquals(0.5, ratios1.get(0).gcContent(), 0.0001);
        assertEquals(72001, ratios1.get(1).Position);
        assertEquals(10.0, ratios1.get(1).ratio(), 0.0001);
        assertEquals(0.5, ratios1.get(1).gcContent(), 0.0001);
        assertEquals(122001, ratios1.get(2).Position);
        assertEquals(172001, ratios1.get(3).Position);

        assertEquals(roundTripped, result);
    }

    @Test
    public void consolidatedRegionsAreLimitedInExtent()
    {
        ListMultimap<Chromosome, BamRatio> ratios = ArrayListMultimap.create();
        for (int i = 0; i < 200; i++)
        {
            int position = i * 100_000 + 1;
            if (i % 10 == 0)
            {
                ratios.put(_1, br(_1, position, 10.0, 0.5, true));
            }
            else
            {
                ratios.put(_1, br(_1, position, -1.0, 0.4, true));
            }
        }

        int consolidationCount = ResultsConsolidator.calcConsolidationCount(8.0);
        LowCoverageConsolidator consolidator = new LowCoverageConsolidator(consolidationCount);
        ListMultimap<Chromosome, BamRatio> result = consolidator.consolidate(ratios);

        Table asTable = toTable(ratios);

        Multimap<String, LowCovBucket> consolidated = LowCoverageRatioMapper.calcConsolidateBuckets(asTable, 8.0);
        assert consolidated != null;
        Table consolidatedTable = new LowCoverageRatioMapper(consolidated, Codec).mapRatios(asTable);
        ListMultimap<Chromosome, BamRatio> roundTripped = toRatios(consolidatedTable);

        assertEquals(1, roundTripped.keySet().size());
        List<BamRatio> ratios1 = roundTripped.get(_1);
        assertEquals(7, ratios1.size());
        assertEquals(1_000_001, ratios1.get(0).Position);
        assertEquals(4_000_001, ratios1.get(1).Position);
        assertEquals(7_000_001, ratios1.get(2).Position);
        assertEquals(10_000_001, ratios1.get(3).Position);
        assertEquals(13_000_001, ratios1.get(4).Position);
        assertEquals(16_000_001, ratios1.get(5).Position);
        assertEquals(18_500_001, ratios1.get(6).Position);

        assertEquals(roundTripped, result);
    }

    private Table toTable(ListMultimap<Chromosome, BamRatio> ratios)
    {
        Table table = Table.create(
                StringColumn.create(CobaltColumns.CHROMOSOME),
                LongColumn.create(ENCODED_CHROMOSOME_POS),
                DoubleColumn.create(RATIO),
                DoubleColumn.create(READ_GC_CONTENT),
                BooleanColumn.create("isAutosome"),
                IntColumn.create(POSITION));

        ratios.keySet().stream().sorted().forEach(chromosome -> {
            List<BamRatio> ratiosForChromosome = ratios.get(chromosome);
            ratiosForChromosome.forEach(bamRatio ->
            {
                Row bucketRow = table.appendRow();
                bucketRow.setString(CobaltColumns.CHROMOSOME, chromosome.contig());
                bucketRow.setLong(ENCODED_CHROMOSOME_POS, Codec.encodeChromosomePosition(chromosome.contig(), bamRatio.Position));
                bucketRow.setDouble(RATIO, bamRatio.ratio());
                bucketRow.setDouble(READ_GC_CONTENT, bamRatio.gcContent());
                bucketRow.setInt(POSITION, bamRatio.Position);
                bucketRow.setBoolean("isAutosome", true);
            });
        });
        return table;
    }

    private ListMultimap<Chromosome, BamRatio>  toRatios(Table table)
    {
        ArrayListMultimap<Chromosome, BamRatio> ratios = ArrayListMultimap.create();
        table.forEach(row -> {
            Chromosome chr = HumanChromosome.fromString(row.getString(CobaltColumns.CHROMOSOME));
            DepthReading depthReading = new DepthReading(chr.contig(), row.getInt(POSITION), row.getDouble(RATIO), row.getDouble(READ_GC_CONTENT) );
            BamRatio ratio = new BamRatio(chr, depthReading, true);
            ratios.put(chr, ratio);
        });
        return ratios;
    }
}
