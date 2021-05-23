package com.hartwig.hmftools.amber;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Random;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.ModifiableBaseDepth;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

public class BaseDepthIntersectFilterTest
{
    private BaseDepth first;
    private BaseDepth second;
    private BaseDepth third;
    private BaseDepthIntersectFilter victim;
    private ListMultimap<Chromosome, BaseDepth> depth;

    @Before
    public void setup()
    {
        final Random random = new Random();
        first = createRandom("1", random);
        second = createRandom("2", random);
        third = createRandom("3", random);
        victim = new BaseDepthIntersectFilter();

        depth = ArrayListMultimap.create();
        depth.put(HumanChromosome._1, first);
        depth.put(HumanChromosome._2, second);
        depth.put(HumanChromosome._3, third);
    }

    @Test
    public void testSingleSample()
    {
        victim.additional(depth.values());
        assertTrue(victim.test(first));
        assertTrue(victim.test(second));
        assertTrue(victim.test(third));
    }

    @Test
    public void testTwoSamples()
    {
        victim.additional(Lists.newArrayList(first, second));
        victim.additional(Lists.newArrayList(third, second));

        assertFalse(victim.test(first));
        assertTrue(victim.test(second));
        assertFalse(victim.test(third));
    }

    @Test
    public void testNoIntersect()
    {
        victim.additional(Lists.newArrayList(first, second));
        victim.additional(Lists.newArrayList(third, second));
        victim.additional(Lists.newArrayList(first, third));

        assertFalse(victim.test(first));
        assertFalse(victim.test(second));
        assertFalse(victim.test(third));
    }

    private BaseDepth createRandom(@NotNull final String chromosome, @NotNull final Random random)
    {
        return ModifiableBaseDepth.create()
                .setChromosome(chromosome)
                .setPosition(random.nextLong())
                .setRef(BaseDepth.Base.A)
                .setAlt(BaseDepth.Base.T)
                .setReadDepth(random.nextInt())
                .setRefSupport(random.nextInt())
                .setAltSupport(random.nextInt())
                .setIndelCount(0);
    }

}
