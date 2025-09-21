package com.hartwig.hmftools.cobalt.targeted;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.count.ReadDepth;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.junit.Assert;
import org.junit.Test;
import org.mockito.Mockito;

public class TargetedRegionEnricherTest
{
    TargetRegionEnricher enricher;
    @Test
    public void constructorTest()
    {
        ListMultimap<Chromosome, TargetRegionEnrichment> initialData = ArrayListMultimap.create();
        initialData.put(_1, tre(_1, 5_001, 0.5));
        initialData.put(_1, tre(_1, 6_001, 0.6));
        initialData.put(_1, tre(_1, 7_001, 0.5));
        initialData.put(_1, tre(_1, 10_001, 0.5));
        initialData.put(_1, tre(_1, 11_001, 0.4));
        initialData.put(_2, tre(_2, 3_001, 0.52));
        initialData.put(_2, tre(_2, 4_001, 0.62));
        TargetRegionEnricher.ChromosomeData chromosomeData = Mockito.mock();
        Mockito.when(chromosomeData.length(Mockito.eq(_1))).thenReturn(15_000);
        Mockito.when(chromosomeData.length(Mockito.eq(_2))).thenReturn(10_000);
        Mockito.when(chromosomeData.length(Mockito.eq(_3))).thenReturn(30_000);

        enricher = new TargetRegionEnricher(initialData, chromosomeData);
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

    void check(double expected, Chromosome chromosome, int position)
    {
        assertEquals(expected, enricher.enrichmentQuotient(chromosome, rd(chromosome, position)), 0.0001);
    }

    private ReadDepth rd(Chromosome chromosome, int position)
    {
        return new ReadDepth(chromosome.contig(), position, 100, 0.5);
    }

    private TargetRegionEnrichment tre(Chromosome chromosome, int position, double enrichment)
    {
        return new TargetRegionEnrichment(chromosome, position, enrichment);
    }
}
