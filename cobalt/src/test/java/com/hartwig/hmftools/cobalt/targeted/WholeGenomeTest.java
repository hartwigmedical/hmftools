package com.hartwig.hmftools.cobalt.targeted;

import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._1;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._2;
import static com.hartwig.hmftools.common.genome.chromosome.HumanChromosome._3;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.cobalt.calculations.DoNothingNormaliser;
import com.hartwig.hmftools.cobalt.count.DepthReading;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;

import org.junit.Before;
import org.junit.Test;

public class WholeGenomeTest
{
    WholeGenome scope;

    @Before
    public void setup(){
        scope = new WholeGenome();
    }

    @Test
    public void enrichmentQuotientTest()
    {
        check(_1, 1);
        check(_1, 1001);
        check(_1, 10_001);
        check(_1, 11_001);
        check(_2, 3001);
        check(_2, 4001);
        check(_3, 1);
        check(_3, 9001);
    }

    @Test
    public void onTargetTest()
    {
        assertTrue(scope.onTarget(_1, 4_999));
        assertTrue(scope.onTarget(_1, 5_000));
        assertTrue(scope.onTarget(_1, 8_000));
        assertTrue(scope.onTarget(_3, 5_001));
    }

    @Test
    public void finalNormaliserTest()
    {
        assertTrue(scope.finalNormaliser() instanceof DoNothingNormaliser);
    }

    void check(Chromosome chromosome, int position)
    {
        assertEquals(1.0, scope.enrichmentQuotient(chromosome, rd(chromosome, position)), 0.0001);
    }

    private DepthReading rd(Chromosome chromosome, int position)
    {
        return new DepthReading(chromosome.contig(), position, 100, 0.5);
    }
}
