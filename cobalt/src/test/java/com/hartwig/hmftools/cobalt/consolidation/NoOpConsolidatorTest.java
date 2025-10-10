package com.hartwig.hmftools.cobalt.consolidation;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.calculations.BamRatio;
import com.hartwig.hmftools.cobalt.calculations.CalculationsTestBase;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.junit.Test;

public class NoOpConsolidatorTest extends CalculationsTestBase
{
    @Test
    public void consolidateTest()
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
        NoOpConsolidator consolidator = new NoOpConsolidator();
        ListMultimap<Chromosome, BamRatio> result = consolidator.consolidate(ratios);

        assertEquals(result, ratios);
    }
}
