package com.hartwig.hmftools.cobalt.targeted;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.calculations.NoOpReadDepthStatisticsNormaliser;
import com.hartwig.hmftools.cobalt.calculations.UnityNormaliser;
import com.hartwig.hmftools.cobalt.consolidation.NoOpConsolidator;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.junit.Before;
import org.junit.Test;
import org.mockito.Mockito;

public class TargetRegionsTest
{
    TargetRegions enricher;

    @Before
    public void setup()
    {
        ListMultimap<Chromosome, TargetRegionEnrichment> initialData = ArrayListMultimap.create();
        initialData.put(_1, tre(_1, 5_001, 0.5));
        initialData.put(_1, tre(_1, 6_001, 0.6));
        initialData.put(_1, tre(_1, 7_001, 0.5));
        initialData.put(_1, tre(_1, 10_001, 0.5));
        initialData.put(_1, tre(_1, 11_001, 0.4));
        initialData.put(_2, tre(_2, 3_001, 0.52));
        initialData.put(_2, tre(_2, 4_001, 0.62));
        TargetRegions.ChromosomeData chromosomeData = Mockito.mock();
        Mockito.when(chromosomeData.length(Mockito.eq(_1))).thenReturn(15_000);
        Mockito.when(chromosomeData.length(Mockito.eq(_2))).thenReturn(10_000);
        Mockito.when(chromosomeData.length(Mockito.eq(_3))).thenReturn(30_000);

        enricher = new TargetRegions(initialData, chromosomeData);
    }

    @Test
    public void enrichmentQuotientTest()
    {
        check(-1.0, _1, 1);
        check(-1.0, _1, 1001);
        check(-1.0, _1, 2001);
        check(-1.0, _1, 3001);
        check(-1.0, _1, 4001);
        check(0.5, _1, 5001);
        check(0.6, _1, 6001);
        check(0.5, _1, 7001);
        check(-1.0, _1, 8001);
        check(-1.0, _1, 9001);
        check(0.5, _1, 10_001);
        check(0.4, _1, 11_001);
        check(-1.0, _1, 12001);
        check(-1.0, _1, 13001);
        check(-1.0, _1, 14001);
        check(-1.0, _2, 1);
        check(-1.0, _2, 1001);
        check(-1.0, _2, 2001);
        check(0.52, _2, 3001);
        check(0.62, _2, 4001);
        check(-1.0, _2, 5001);
        check(-1.0, _2, 9001);

        check(-1.0, _3, 1);
        check(-1.0, _3, 9001);
    }

    @Test
    public void onTargetTest()
    {
        assertFalse(enricher.onTarget(_1, 4_999));
        assertTrue(enricher.onTarget(_1, 5_000));
        assertTrue(enricher.onTarget(_1, 5_001));
        assertTrue(enricher.onTarget(_1, 5_031));
        assertTrue(enricher.onTarget(_1, 5_999));
        assertTrue(enricher.onTarget(_1, 7_999));
        assertFalse(enricher.onTarget(_1, 8_000));

        assertFalse(enricher.onTarget(_3, 5_001));
    }

    @Test
    public void resultsConsolidator()
    {
        assertTrue(enricher.resultsConsolidator(0.67) instanceof NoOpConsolidator);
        assertTrue(enricher.resultsConsolidator(670.0) instanceof NoOpConsolidator);
    }

    @Test
    public void finalNormaliserTest()
    {
        assertTrue(enricher.finalNormaliser() instanceof UnityNormaliser);
    }

    @Test
    public void medianByMeanNormaliserTest()
    {
        assertTrue(enricher.medianByMeanNormaliser() instanceof NoOpReadDepthStatisticsNormaliser);
    }

    void check(double expected, Chromosome chromosome, int position)
    {
        assertEquals(expected, enricher.enrichmentQuotient(chromosome, rd(chromosome, position)), 0.0001);
    }

    private DepthReading rd(Chromosome chromosome, int position)
    {
        return new DepthReading(chromosome.contig(), position, 100, 0.5);
    }

    private TargetRegionEnrichment tre(Chromosome chromosome, int position, double enrichment)
    {
        return new TargetRegionEnrichment(chromosome, position, enrichment);
    }
}
