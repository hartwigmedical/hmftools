package com.hartwig.hmftools.common.amber;

import static org.junit.Assert.assertEquals;

import java.util.Random;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Before;
import org.junit.Test;

public class BaseDepthIntersectTest {

    private BaseDepth first;
    private BaseDepth second;
    private BaseDepth third;
    private BaseDepthIntersect victim;
    private ListMultimap<Chromosome, BaseDepth> depth;

    @Before
    public void setup() {
        final Random random = new Random();
        first = TumorContaminationFileTest.createRandom("1", random);
        second = TumorContaminationFileTest.createRandom("2", random);
        third = TumorContaminationFileTest.createRandom("3", random);
        victim = new BaseDepthIntersect();

        depth = ArrayListMultimap.create();
        depth.put(HumanChromosome._1, first);
        depth.put(HumanChromosome._2, second);
        depth.put(HumanChromosome._3, third);

    }

    @Test
    public void testSingleSample() {
        victim.add(depth.values());
        ListMultimap<Chromosome, BaseDepth> result = victim.intersect(depth);
        assertEquals(3, result.size());
    }

    @Test
    public void testTwoSamples() {
        victim.add(Lists.newArrayList(first, second));
        victim.add(Lists.newArrayList(third, second));

        ListMultimap<Chromosome, BaseDepth> result = victim.intersect(depth);
        assertEquals(1, result.size());
        assertEquals(second, result.get(HumanChromosome._2).get(0));
    }

    @Test
    public void testNoIntersect() {
        victim.add(Lists.newArrayList(first, second));
        victim.add(Lists.newArrayList(third, second));
        victim.add(Lists.newArrayList(first, third));

        ListMultimap<Chromosome, BaseDepth> result = victim.intersect(depth);
        assertEquals(0, result.size());
    }

}
