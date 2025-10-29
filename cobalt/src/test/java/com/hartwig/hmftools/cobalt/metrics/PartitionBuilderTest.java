package com.hartwig.hmftools.cobalt.metrics;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.ListMultimap;
import com.hartwig.hmftools.cobalt.ChromosomeData;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;

import org.junit.Test;

public class PartitionBuilderTest
{
    @Test
    public void partitionChromosomesTest()
    {
        List<ChromosomeData> chromosomes = new ArrayList<>();
        chromosomes.add(new ChromosomeData("1", 249250621));
        chromosomes.add(new ChromosomeData("3", 198022430));
        chromosomes.add(new ChromosomeData("Y", 59_373_566));

        ListMultimap<Chromosome, Partition> result = PartitionBuilder.partitionChromosomes(chromosomes, 10_000_000, 1_000_000);
        assertEquals(3, result.keySet().size());

        List<Partition> partitionsY = result.get(HumanChromosome._Y);
        assertEquals(5, partitionsY.size());
        Partition partitionY0 = partitionsY.get(0);
        assertEquals(HumanChromosome._Y.shortName(), partitionY0.chromosome());
        assertEquals(1, partitionY0.start());
        assertEquals(10_000_000, partitionY0.end());
        assertEquals(10, partitionY0.TargetRegions.size());
        assertEquals(1, partitionY0.TargetRegions.get(0).start());
        assertEquals(1_000_000, partitionY0.TargetRegions.get(0).end());
    }

    @Test
    public void createPartitionsTest()
    {
        PartitionBuilder builder = new PartitionBuilder(1_000_000, 10_000);
        ChromosomeData chromosomeData = new ChromosomeData("21", 48_129_895);
        List<Partition> partitions = builder.createPartitions(chromosomeData);
        assertEquals(48, partitions.size());

        assertEquals(1, partitions.get(0).start());
        assertEquals(1_000_000, partitions.get(0).end());
        assertEquals("21", partitions.get(0).Chromosome);

        assertEquals(10_000_001, partitions.get(10).start());
        assertEquals(11_000_000, partitions.get(10).end());
        assertEquals("21", partitions.get(10).Chromosome);

        assertEquals(47_000_001, partitions.get(47).start());
        assertEquals(48_000_000, partitions.get(47).end());
        assertEquals("21", partitions.get(47).Chromosome);
    }
}
