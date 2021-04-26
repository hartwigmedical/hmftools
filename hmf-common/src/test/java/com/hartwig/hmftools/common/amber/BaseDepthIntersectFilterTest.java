package com.hartwig.hmftools.common.amber;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Random;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Before;
import org.junit.Test;

public class BaseDepthIntersectFilterTest {

    private BaseDepth first;
    private BaseDepth second;
    private BaseDepth third;
    private BaseDepthIntersectFilter victim;
    private ListMultimap<Chromosome, BaseDepth> depth;

    @Before
    public void setup() {
        final Random random = new Random();
        first = TumorContaminationFileTest.createRandom("1", random);
        second = TumorContaminationFileTest.createRandom("2", random);
        third = TumorContaminationFileTest.createRandom("3", random);
        victim = new BaseDepthIntersectFilter();

        depth = ArrayListMultimap.create();
        depth.put(HumanChromosome._1, first);
        depth.put(HumanChromosome._2, second);
        depth.put(HumanChromosome._3, third);
    }

    @Test
    public void testSingleSample() {
        victim.additional(depth.values());
        assertTrue(victim.test(first));
        assertTrue(victim.test(second));
        assertTrue(victim.test(third));
    }

    @Test
    public void testTwoSamples() {
        victim.additional(Lists.newArrayList(first, second));
        victim.additional(Lists.newArrayList(third, second));

        assertFalse(victim.test(first));
        assertTrue(victim.test(second));
        assertFalse(victim.test(third));
    }

    @Test
    public void testNoIntersect() {
        victim.additional(Lists.newArrayList(first, second));
        victim.additional(Lists.newArrayList(third, second));
        victim.additional(Lists.newArrayList(first, third));

        assertFalse(victim.test(first));
        assertFalse(victim.test(second));
        assertFalse(victim.test(third));
    }
}
