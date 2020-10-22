package com.hartwig.hmftools.common.cobalt;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Test;

public class MedianRatioFactoryTest {
    @Test
    public void testHg19() {
        ReadRatio hg19Ratio = ImmutableReadRatio.builder().chromosome("1").position(1).ratio(1).build();
        Multimap<Chromosome, ReadRatio> map = ArrayListMultimap.create();
        map.put(HumanChromosome._1, hg19Ratio);
        List<MedianRatio> victim = MedianRatioFactory.createFromReadRatio(map);
        assertEquals(1, victim.size());
        assertEquals("1", victim.get(0).chromosome());
    }

    @Test
    public void testHg38() {
        ReadRatio hg19Ratio = ImmutableReadRatio.builder().chromosome("chr1").position(1).ratio(1).build();
        Multimap<Chromosome, ReadRatio> map = ArrayListMultimap.create();
        map.put(HumanChromosome._1, hg19Ratio);
        List<MedianRatio> victim = MedianRatioFactory.createFromReadRatio(map);
        assertEquals(1, victim.size());
        assertEquals("chr1", victim.get(0).chromosome());
    }

}
