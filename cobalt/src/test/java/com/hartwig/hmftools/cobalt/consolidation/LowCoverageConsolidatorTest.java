package com.hartwig.hmftools.cobalt.consolidation;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.calculations.BamRatio;
import com.hartwig.hmftools.cobalt.calculations.CalculationsTestBase;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.junit.Test;

public class LowCoverageConsolidatorTest extends CalculationsTestBase
{
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

        assertEquals(2, result.keySet().size());
        List<BamRatio> ratios1 = result.get(_1);
        List<BamRatio> ratios2 = result.get(_2);
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
    }

    @Test
    public void useConsolidationCount()
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
        int consolidationCount = ResultsConsolidator.calcConsolidationCount(4.0);
        LowCoverageConsolidator consolidator = new LowCoverageConsolidator(consolidationCount);
        ListMultimap<Chromosome, BamRatio> result = consolidator.consolidate(ratios);

        assertEquals(2, result.keySet().size());
        List<BamRatio> ratios1 = result.get(_1);
        List<BamRatio> ratios2 = result.get(_2);
        assertEquals(2, ratios1.size());
        assertEquals(2, ratios2.size());
        assertEquals(9001, ratios1.get(0).Position);
        assertEquals(10.095, ratios1.get(0).ratio(), 0.0001);
        assertEquals(0.495, ratios1.get(0).gcContent(), 0.0001);
        assertEquals(9001, ratios2.get(0).Position);
        assertEquals(10.095, ratios2.get(0).ratio(), 0.0001);
        assertEquals(0.495, ratios2.get(0).gcContent(), 0.0001);
        assertEquals(30001, ratios1.get(1).Position);
        assertEquals(10.295, ratios1.get(1).ratio(), 0.0001);
        assertEquals(0.695, ratios1.get(1).gcContent(), 0.0001);
        assertEquals(30001, ratios2.get(1).Position);
        assertEquals(10.295, ratios2.get(1).ratio(), 0.0001);
        assertEquals(0.695, ratios2.get(1).gcContent(), 0.0001);
    }

    @Test
    public void maskedRatiosAreSkippedInConsolidation()
    {
        ListMultimap<Chromosome, BamRatio> ratios = ArrayListMultimap.create();
        for(int i = 0; i < 200; i++)
        {
            int position = i * 1000 + 1;
            if(i % 5 == 0)
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

        assertEquals(1, result.keySet().size());
        List<BamRatio> ratios1 = result.get(_1);
        assertEquals(4, ratios1.size());
        assertEquals(23001, ratios1.get(0).Position);
        assertEquals(10.0, ratios1.get(0).ratio(), 0.0001);
        assertEquals(0.5, ratios1.get(0).gcContent(), 0.0001);
        assertEquals(72001, ratios1.get(1).Position);
        assertEquals(10.0, ratios1.get(1).ratio(), 0.0001);
        assertEquals(0.5, ratios1.get(1).gcContent(), 0.0001);
        assertEquals(122001, ratios1.get(2).Position);
        assertEquals(172001, ratios1.get(3).Position);
    }

    @Test
    public void consolidatedRegionsAreLimitedInExtent()
    {
        ListMultimap<Chromosome, BamRatio> ratios = ArrayListMultimap.create();
        for(int i = 0; i < 200; i++)
        {
            int position = i * 100_000 + 1;
            if(i % 10 == 0)
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

        assertEquals(1, result.keySet().size());
        List<BamRatio> ratios1 = result.get(_1);
        assertEquals(7, ratios1.size());
        assertEquals(1_000_001, ratios1.get(0).Position);
        assertEquals(4_000_001, ratios1.get(1).Position);
        assertEquals(7_000_001, ratios1.get(2).Position);
        assertEquals(10_000_001, ratios1.get(3).Position);
        assertEquals(13_000_001, ratios1.get(4).Position);
        assertEquals(16_000_001, ratios1.get(5).Position);
        assertEquals(18_500_001, ratios1.get(6).Position);
    }

    @Test
    public void boundariesAreReused()
    {
        ListMultimap<Chromosome, BamRatio> tumorRatios = ArrayListMultimap.create();
        for(int i = 0; i < 1000; i++)
        {
            int position = i * 1000 + 1;
            if(i % 10 == 0)
            {
                tumorRatios.put(_1, br(_1, position, 10.0, 0.5, true));
                tumorRatios.put(_2, br(_2, position, 12.0, 0.52, true));
            }
            else
            {
                tumorRatios.put(_1, br(_1, position, -1.0, 0.4, true));
                tumorRatios.put(_2, br(_2, position, -1.0, 0.42, true));
            }
        }

        int consolidationCount = ResultsConsolidator.calcConsolidationCount(8.0);
        LowCoverageConsolidator consolidator = new LowCoverageConsolidator(consolidationCount);
        ListMultimap<Chromosome, BamRatio> tumCons = consolidator.consolidate(tumorRatios);

        ListMultimap<Chromosome, BamRatio> referenceRatios = ArrayListMultimap.create();
        for(int i = 0; i < 1000; i++)
        {
            int position = i * 1000 + 1;
            referenceRatios.put(_1, br(_1, position, 11.0, 0.51, true));
            referenceRatios.put(_2, br(_2, position, 13.0, 0.53, true));
        }
        ListMultimap<Chromosome, BamRatio> refCons = consolidator.consolidate(referenceRatios);

        List<Integer> consolidatedTumorRatioPositionsChr1 = tumCons.get(_1).stream().map(bamRatio -> bamRatio.Position).toList();
        List<Integer> consolidatedTumorRatioPositionsChr2 = tumCons.get(_2).stream().map(bamRatio -> bamRatio.Position).toList();
        List<Integer> consolidatedRefRatioPositionsChr1 = refCons.get(_1).stream().map(bamRatio -> bamRatio.Position).toList();
        List<Integer> consolidatedRefRatioPositionsChr2 = refCons.get(_2).stream().map(bamRatio -> bamRatio.Position).toList();
        assertEquals(consolidatedTumorRatioPositionsChr1, consolidatedRefRatioPositionsChr1);
        assertEquals(consolidatedTumorRatioPositionsChr2, consolidatedRefRatioPositionsChr2);
    }
}
