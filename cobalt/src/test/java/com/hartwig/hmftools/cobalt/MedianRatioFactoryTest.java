package com.hartwig.hmftools.cobalt;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.cobalt.MedianRatioFactory;
import com.hartwig.hmftools.common.cobalt.ImmutableReadRatio;
import com.hartwig.hmftools.common.cobalt.MedianRatio;
import com.hartwig.hmftools.common.cobalt.ReadRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Test;

public class MedianRatioFactoryTest {

    @Test
    public void testV37() {
        ReadRatio v37Ratio = ImmutableReadRatio.builder().chromosome("1").position(1).ratio(1).build();
        Multimap<Chromosome, ReadRatio> map = ArrayListMultimap.create();
        map.put(HumanChromosome._1, v37Ratio);
        List<MedianRatio> victim = MedianRatioFactory.createFromReadRatio(map);
        assertEquals(1, victim.size());
        assertEquals("1", victim.get(0).chromosome());
    }

    @Test
    public void testV38() {
        ReadRatio v38Ratio = ImmutableReadRatio.builder().chromosome("chr1").position(1).ratio(1).build();
        Multimap<Chromosome, ReadRatio> map = ArrayListMultimap.create();
        map.put(HumanChromosome._1, v38Ratio);
        List<MedianRatio> victim = MedianRatioFactory.createFromReadRatio(map);
        assertEquals(1, victim.size());
        assertEquals("chr1", victim.get(0).chromosome());
    }
}
